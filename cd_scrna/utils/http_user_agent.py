# Helper imported by cd_scrna/14_treg_pathway_tf_activity.py (section F4G_pathways).
"""Send a browser User-Agent with requests.get so that decoupler can download its resources.

Zenodo answers the default python-requests User-Agent with HTTP 403. Importing this module wraps requests.get (and the
reference used by decoupler's download module) so that every call carries a browser User-Agent; decoupler's own
dc.op.collectri / progeny / hallmark loaders then work unchanged."""
import requests, functools
_orig_get = requests.get
@functools.wraps(_orig_get)
def _get(url, *a, **k):
    h = dict(k.pop('headers', {}) or {}); h.setdefault('User-Agent', 'Mozilla/5.0 (X11; Linux x86_64)')
    return _orig_get(url, *a, headers=h, **k)
requests.get = _get
import decoupler._download as _d
_d.requests.get = _get
