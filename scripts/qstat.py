#!/usr/bin/env python
import sys
if not ((2, 7) <= sys.version_info[:2]):
    sys.exit('Python 2.7 and 3 ares supported '
             '(you are running %d.%d.%d)' %
             (sys.version_info[0], sys.version_info[1], sys.version_info[2]))

from itertools import repeat, count
import subprocess
import sys
import getpass


example_lines = '''job-ID     prior   name       user         state submit/start at     queue                          jclass                         slots ja-task-ID
------------------------------------------------------------------------------------------------------------------------------------------------
      9606 0.86405 VFS_BxAUAg klpf990      r     03/18/2015 11:30:23 batch.q@bn0204                                                    5
       Full jobname:     VFS_BxAUAg=_bcbio
       Master Queue:     batch.q@bn0204
       Requested PE:     smp 5
       Granted PE:       smp 5
       Hard Resources:
       Soft Resources:
       Hard requested queues: batch.q
       Predecessor Jobs (request): _
       Binding:          NONE
      9670 0.00000 TARGQC_0G6 klpf990      hqw   03/18/2015 11:37:00                                                                   4
       Full jobname:     TARGQC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup
       Requested PE:     smp 4
       Hard Resources:
       Soft Resources:
       Hard requested queues: batch.q
       Predecessor Jobs (request): TC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Acrometrix-ready, NC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Acrometrix-ready, QM_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Acrometrix-ready, TC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Colo829-ready, NC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Colo829-ready, QM_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Colo829-ready, TC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Horizon-ready, NC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Horizon-ready, QM_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Horizon-ready, TC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Promega-ready, NC_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Promega-ready, QM_0G6ozA=_IDT-PanCancer_AZ1-IDT-Evaluation_nodup_Promega-ready
       Predecessor Jobs: 9662, 9668, 9244, 9544
       Binding:          NONE'''.split('\n')


rows = []
header_fields = []
cur_fields, full_name, pred_jobs = None, '', ''

def add_row(cur_fields, full_name, pred_jobs):
    if len(cur_fields) >= 9:
        job_id, prior, name, user, state, date, time, queue, slots = cur_fields[:9]
        rows.append([job_id, prior, full_name, user, state, date, time, queue, slots, pred_jobs])
    else:
        if len(cur_fields) >= 8:
            job_id, prior, name, user, state, date, time, slots = cur_fields[:8]
            queue = ''
            rows.append([job_id, prior, full_name, user, state, date, time, queue, slots, pred_jobs])
        else:
            sys.stderr.write('Line with small number of fields: ' + str(cur_fields) + '\n')
    return None, '', ''


un = getpass.getuser()
cmdl = ['qstat', '-r']
if len(sys.argv) > 1 and sys.argv[1] == '-a':
    un = ''
    cmdl.extend(['-u', '"*"'])

f = subprocess.Popen(cmdl, stdout=subprocess.PIPE).stdout
for i, l in enumerate(f):
    if i == 0:
        l = l.replace('submit/start', 'submit/date')
        header_fields = [f for f in l.split() if f not in ['jclass', 'ja-task-ID']]
        rows.append(header_fields + ['pred_jobs'])  # job-ID  prior  name  user  state  submit/start  at  queue  slots
        continue
    if i == 1:  # ----------------...
        continue

    l = l.strip()
    if un in l.strip():  # 9606 0.86405 VFS_BxAUAg klpf990      r     03/18/2015 11:30:23 batch.q@bn0204
        if cur_fields:
            cur_fields, full_name, pred_jobs = add_row(cur_fields, full_name, pred_jobs)
        cur_fields = l.split()

    elif cur_fields:
        if l.startswith('Full jobname: '):
            full_name = l.split(': ')[1].strip()

        elif l.startswith('Predecessor Jobs: '):
            pred_jobs = ' '.join(l.split(': ')[1].strip().split(', '))

if cur_fields:
    add_row(cur_fields, full_name, pred_jobs)


if len(rows) > 0:
    col_widths = repeat(0)
    try:
        from itertools import izip as zip
    except:
        pass
    for row in rows:
        col_widths = [max(len(v), w) for v, w in zip(row[:-1], col_widths)]

    line = ''
    for i, row in enumerate(rows):
        if i == 1:
            sys.stdout.write('-' * len(line) + '\n')
        for val, w in zip(row[:-1], col_widths):
            cell = val + (' ' * (w - len(val) + 2))
            sys.stdout.write(cell)
            line += cell
        sys.stdout.write(row[-1])
        line += row[-1]
        sys.stdout.write('\n')