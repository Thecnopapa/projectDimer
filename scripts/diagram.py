
import sys, os

from pyflowchart import Flowchart, output_html


print(sys.argv)

if len(sys.argv)>2:
    in_file = sys.argv [2]

    with open('scripts/{}.py'.format(in_file)) as f:
        code = f.read()

    fc = Flowchart.from_code(code, inner=True, simplify=True, conds_align=False)
    print(fc.flowchart())
    os.makedirs('diagrams', exist_ok=True)
    output_html('diagrams/{}.html'.format(in_file), "", flowchart=fc.flowchart())
    print("Diagram saved to:", 'diagrams/{}.html'.format(in_file))

    # output flowchart code.("userdb")
else:
    print("No input file")







