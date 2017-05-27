#!/usr/bin/env python

import os
import sys
import base64
import mimetypes
import traceback
from bs4 import BeautifulSoup


def execute(html_fpath):
    if not os.path.isfile(html_fpath):
        return None
    b64_html_fpath = os.path.splitext(html_fpath)[0] + '.b64html'

    data = _html_to_dataurl(html_fpath)

    with open(b64_html_fpath, 'wt') as fp:
        fp.write(data)
    print('Saved to ' + b64_html_fpath)


def _html_to_dataurl(html_fpath):
    comp_html_unicode = _compact_html(html_fpath)
    return _convert_text(comp_html_unicode.encode('utf-8'), mime='text/html')


def _compact_html(html_fpath):
    """ 1. Embed images in base64 format
        2. Replace anchors with javascript
    """
    with open(html_fpath) as f:
        t = f.read()
        try:
            html_unicode = unicode(t, 'utf-8')
        except TypeError:
            html_unicode = unicode(t)
        except NameError:  # py3
            html_unicode = t
    base_dir = os.path.split(html_fpath)[0]
    soup = BeautifulSoup(html_unicode, 'lxml')
    for img in soup.findAll('img'):
        img_src = img['src']
        durl_img = _convert_file_contents(os.path.join(base_dir, img['src']))
        img['src'] = durl_img

    js = "javascript: void(0); document.getElementById('%s').scrollIntoView(true);"
    for anchor in soup.findAll('a'):
        if 'href' not in anchor:
            continue
        if anchor['href'].startswith('#'):
            anchor['href'] = js % anchor['href'][1:]
        else:
            del anchor['href']
    return soup.prettify()


def _convert_file_contents(fpath, mime=None):
    print('Encoding file ' + fpath)
    if not os.path.isfile(fpath):
        sys.stderr.write('File ' + fpath + ' does not exists\n')
        return ''
    if not mime:
        mimetypes.init()
        mime, enc = mimetypes.guess_type(os.path.join('file://', fpath))
        if mime is None:
            raise Exception('rfc2397: failed to determine file type')
        print('Mime: ' + str(mime) + ', encoding name: ' + str(enc))

    with open(fpath, 'rb') as fp:
        data = fp.read()
    b64_string = base64.b64encode(data)
    return 'data:%s;base64,%s' % (mime, b64_string)


def _convert_text(text, mime):
    print('Encoding ' + 'Mime: ' + str(mime))
    b64_string = base64.b64encode(text)
    return 'data:%s;base64,%s' % (mime, b64_string)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.stderr.write('Please, provide input HTML file. Usage: ' + __file__ + ' file.html\n')
        sys.exit(1)
    execute(sys.argv[1])


