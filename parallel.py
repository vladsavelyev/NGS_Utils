import contextlib
import subprocess

from cluster_helper.cluster import ClusterView as CV
from joblib import Parallel, delayed

from Utils.logger import debug
from Utils.utils import is_local


class ParallelCfg:
    def __init__(self,
                 scheduler=None,
                 queue=None,
                 resources=None,
                 threads=None,
                 tag='targqc'):
        self.scheduler = scheduler
        self.queue = queue
        self.threads = threads or 1
        self.extra_params = dict()
    
        self.extra_params['run_local'] = is_local()
        if tag:
            self.extra_params['tag'] = tag
        if resources:
            self.extra_params['resources'] = resources

    def num_jobs(self, n_samples):
        return min(self.threads, n_samples)

    def cores_per_job(self, n_jobs):
        return max(1, self.threads / n_jobs)

    def get_cluster_params(self, n_samples):
        return dict(
            scheduler=self.scheduler,
            queue=self.queue,
            num_jobs=self.num_jobs(n_samples),
            cores_per_job=self.cores_per_job(self.num_jobs(n_samples)),
            extra_params=self.extra_params)


def get_parallel_view(n_samples, parallel_cfg):
    if parallel_cfg.scheduler:
        debug('Starting' + (' test' if is_local() else '') + ' cluster (scheduler: ' + parallel_cfg.scheduler + ', queue: ' + parallel_cfg.queue + ') '
              'using ' + str(parallel_cfg.num_jobs(n_samples)) + ' nodes, ' + str(parallel_cfg.cores_per_job(n_samples)) + ' thread per each sample')
        return ClusterView(n_samples, parallel_cfg)
    else:
        debug('Running locally using ' + str(parallel_cfg.num_jobs(n_samples)) + ' thread(s)')
        return ThreadedView(n_samples, parallel_cfg)


@contextlib.contextmanager
def parallel_view(n_samples, parallel_cfg):
    view = get_parallel_view(n_samples, parallel_cfg)
    try:
        yield view
    finally:
        view.stop()


class BaseView:
    def __init__(self, n_samples, parallel_cfg):
        self.n_samples = n_samples
        self.parallel_cfg = parallel_cfg
        self.num_jobs = parallel_cfg.num_jobs(n_samples)
        self.cores_per_job = parallel_cfg.cores_per_job(n_samples)
        self._view = None

    def run(self, fn, param_lists):
        raise NotImplementedError

    def stop(self):
        raise NotImplementedError


class ClusterView(BaseView):
    def __init__(self, n_samples, parallel_cfg):
        BaseView.__init__(self, n_samples, parallel_cfg)
        self._view = CV(**parallel_cfg.get_cluster_params(n_samples))

    def run(self, fn, param_lists):
        assert self.n_samples == len(param_lists)
        n_params = len(param_lists[0])
        return self._view.view.map(fn, *([params[param_i] for params in param_lists] for param_i in range(n_params)))

    def stop(self):
        self._view.stop()


class ThreadedView(BaseView):
    def __init__(self, n_samples, parallel_cfg):
        BaseView.__init__(self, n_samples, parallel_cfg)
        self._view = Parallel(n_jobs=self.num_jobs)

    def run(self, fn, param_lists):
        assert self.n_samples == len(param_lists)
        return self._view(delayed(fn)(*params) for params in param_lists)

    def stop(self):
        return