
import sys, os

from pyflowchart import Flowchart, output_html
from utilities import *

print(sys.argv)


sprint("Generating diagram...")
if len(sys.argv)>1:
    in_file = sys.argv [1]
    print1("In file:", in_file)
    if len(sys.argv)>2:
        field=sys.argv[2]
    else:
        field=''
    print1("Field:", field)

    with open('scripts/{}.py'.format(in_file)) as f:
        print1('Opening file at: scripts/{}.py'.format(in_file))
        code = f.read()

    fc = Flowchart.from_code(code, field=field, inner=True, simplify=True, conds_align=False)
    #print(fc.flowchart())
    os.makedirs('diagrams', exist_ok=True)

    if field == '':
        field = "all"
    out_folder = "diagrams/{}".format(in_file)
    print(out_folder)
    os.makedirs(out_folder, exist_ok=True)
    out_path = os.path.abspath(os.path.join(out_folder,"{}.html".format(field)))

    output_html(out_path, "", flowchart=fc.flowchart())
    print1("Diagram saved to:", out_path)

    # output flowchart code.("userdb")
else:
    sprint("No input file")







