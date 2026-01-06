import re
import sys


def parse_c_dsm(filepath):
    with open(filepath, "r") as f:
        content = f.read()

    # Extract the dsm_t04_pos initializer
    match = re.search(
        r"const\s+dsm_t\s+dsm_t04_pos\s*=\s*{(.+?)}\s*;", content, re.DOTALL
    )
    if not match:
        print("Could not find dsm_t04_pos in C file")
        return None

    data_str = match.group(1)
    # Extract all integers
    integers = re.findall(r"-?\d+", data_str)
    integers = [int(x) for x in integers]

    if len(integers) != 6 * 6 * 6 * 6:
        print(f"C file parsed {len(integers)} ints, expected {6 * 6 * 6 * 6}")
        return None

    return integers


def generate_rust_table(values):
    print("pub const DSM_T04_POS: DsmTable = [")
    idx = 0
    for q1 in range(6):
        print("    [")
        for q2 in range(6):
            print("        [")
            for t1 in range(6):
                print("            [", end="")
                row_vals = []
                for t2 in range(6):
                    row_vals.append(str(values[idx]))
                    idx += 1
                print(", ".join(row_vals), end="")
                print("],")
            print("        ],")
        print("    ],")
    print("];")


c_vals = parse_c_dsm(
    "/Users/vzr518/Projects/risearch/legacy_c/RIsearch2/src/dsm-t04-guugle.c"
)

if c_vals:
    generate_rust_table(c_vals)
else:
    print("Failed to parse C values", file=sys.stderr)
