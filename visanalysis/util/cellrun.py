# -*- coding: utf-8 -*-
"""Run an analysis script from a shell-script cell and see its output there.

WHY THIS EXISTS. The shell scripts used os.system, which hands the child
process the CONSOLE's stdout. Spyder and the VS Code interactive window only
redirect Python-level sys.stdout, so everything the analysis printed went to
the external console and the cell that launched it looked empty -- the output
had to be read by retyping the command in a terminal.

WHY A MODULE RATHER THAN A DEF IN THE FIRST CELL. A cell-local helper only
exists in a kernel that has run that cell since the file was last edited, so
editing the shell script mid-session makes every later block fail with
`NameError: name 'run' is not defined`. Importing it means a block depends only
on the arguments cell, not on when the environment cell last ran.
"""
import os
import shlex
import subprocess
import sys


def run(script, *args):
    """Run `script` in this interpreter and print its output into this cell.

    sys.executable, not 'python': the child then runs in the SAME interpreter
    as the shell script, rather than whichever python happens to be first on
    PATH -- which may be none at all, and with os.system that failure printed
    to the console and left the cell silent.

    -u keeps the child unbuffered, so lines appear while it works instead of
    arriving in one block at the end; these scripts run for minutes.

    Arguments are passed as a list, so a path containing a space stays one
    argument instead of silently splitting.

    Returns the child's exit code (0 on success).
    """
    cmd = [sys.executable, "-u", str(script)] + [str(a) for a in args]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, bufsize=1)
    for line in p.stdout:
        print(line, end="")
        sys.stdout.flush()
    code = p.wait()
    if code:
        print("!! {} exited with code {}".format(
            os.path.basename(str(script)), code))
    return code


def run_cmd(command):
    """Run a command STRING and print its output into this cell.

    The older shell scripts build one string by concatenation and handed it to
    os.system; this is the drop-in replacement, so converting them is a
    one-token change per call site rather than a rewrite of each command.

    A leading bare `python` is replaced with sys.executable, so the analysis
    runs in the SAME interpreter as the shell script instead of whichever
    python is first on PATH -- which may be none, a failure os.system reported
    only to the console while the cell stayed empty.

    posix=False keeps Windows backslashes intact through the split, but it
    also keeps the QUOTES: `--flag "a b"` splits into the token `"a b"` with
    the quote characters still attached. os.system handed the whole string to
    a shell, which removed them; nothing removes them here, so a quoted path
    reaches the child as `"C:/a b/x.h5"` and fails to open, and a quoted
    number reaches argparse as `"10"` and fails to parse. So strip one
    balanced pair per token -- the single thing a shell did that the
    posix=False split does not. Prefer run() with a list when writing new code.
    """
    parts = shlex.split(str(command), posix=False)
    parts = [p[1:-1] if len(p) > 1 and p[0] == p[-1] and p[0] in '"\'' else p
             for p in parts]
    if parts and os.path.basename(parts[0]).lower().strip('"').startswith("python"):
        parts = parts[1:]
    return run(*parts) if parts else 1
