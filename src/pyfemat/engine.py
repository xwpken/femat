import matlab.engine
from pathlib import Path

_ROOT = Path(__file__).resolve().parent.parent.parent
_eng = None


def _start():
    global _eng
    _eng = matlab.engine.start_matlab()
    _eng.eval(f"addpath(genpath('{_ROOT / 'src'}'));", nargout=0)
    _eng.eval(f"addpath('{_ROOT / 'tests'}');", nargout=0)
    return _eng


def eng():
    if _eng is None:
        return _start()
    return _eng


def put(name, value):
    eng().workspace[name] = value


def execute(cmd):
    eng().eval(cmd, nargout=0)


def get(name):
    return eng().workspace[name]


def call(func, *args, nargout=1):
    args = [float(a) if isinstance(a, int) else a for a in args]
    return eng().feval(func, *args, nargout=nargout)
