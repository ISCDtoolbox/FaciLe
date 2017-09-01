import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__),"../projects/tools"))
print os.path.join(os.path.dirname(os.path.abspath(__file__)),"../tools")
import msh
import executable_paths as exe

def command(cmd, displayOutput=False):
    err = 1
    print "Running the command '" + cmd + "'"

    if displayOutput:
        err = os.system(cmd)
    else:
        err = os.system(cmd + " > tmp_out.txt 2>tmp_err.txt")

    if err:
        print "An error happened while executing:\n"+cmd+"\nLook in tmp_out.txt or tmp_err.txt for info\nExiting..."
        sys.exit()
    else:
        os.system("rm tmp_out.txt tmp_err.txt >/dev/null 2>&1")

if __name__=="__main__":
    command(exe.shell + " ~/source.mesh", displayOutput=True)
