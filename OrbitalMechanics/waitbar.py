import sys

__author__ = 'Ian'


# Found main body of source code @ http://stackoverflow.com/questions/3160699/python-progress-bar

def bar(progress, message):
    barLength = 10  # Modify this to change the length of the progress bar
    status = message
    done = False
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
        done = True
    block = int(round(barLength * progress))
    text = "\rPercent: [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), progress * 100, status)
    sys.stdout.write(text)
    sys.stdout.flush()
    return done
