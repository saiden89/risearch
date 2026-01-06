import re


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
    # Remove newlines and outer braces
    data_str = data_str.replace("\n", "")

    # We want to parse this nested structure. Only numbers matter.
    # The structure is fixed 6x6x6x6.
    # Let's just extract all numbers in order.
    integers = re.findall(r"-?\d+", data_str)
    integers = [int(x) for x in integers]

    if len(integers) != 6 * 6 * 6 * 6:
        print(f"C file parsed {len(integers)} ints, expected {6 * 6 * 6 * 6}")
        return None

    return integers


def parse_rust_dsm(filepath):
    with open(filepath, "r") as f:
        content = f.read()

    match = re.search(
        r"pub\s+const\s+DSM_T04_POS:\s*DsmTable\s*=\s*\[(.+?)\s*\];", content, re.DOTALL
    )
    if not match:
        print("Could not find DSM_T04_POS in Rust file")
        return None

    data_str = match.group(1)  # Capture only the content inside the brackets
    integers = re.findall(r"-?\d+", data_str)
    integers = [int(x) for x in integers]

    if len(integers) != 6 * 6 * 6 * 6:
        print(f"Rust file parsed {len(integers)} ints, expected {6 * 6 * 6 * 6}")
        print(f"First 10: {integers[:10]}")
        print(f"Last 10: {integers[-10:]}")
        # return None # Try to proceed or debug

    return integers


c_vals = parse_c_dsm(
    "/Users/vzr518/Projects/risearch/legacy_c/RIsearch2/src/dsm-t04-guugle.c"
)
rust_vals = parse_rust_dsm("/Users/vzr518/Projects/risearch/src/dsm.rs")

if c_vals and rust_vals:
    print(f"C len: {len(c_vals)}, Rust len: {len(rust_vals)}")

    # Truncate if Rust has 1 extra? Or fix regex?
    # Let's ensure we align them.
    if len(rust_vals) > len(c_vals):
        print("Rust has more values. Checking alignment...")

    length = min(len(c_vals), len(rust_vals))
    diffs = 0
    for i in range(len(c_vals)):
        if c_vals[i] != rust_vals[i]:
            # decode indices
            rem = i
            l = rem % 6
            rem //= 6
            k = rem % 6
            rem //= 6
            j = rem % 6
            rem //= 6
            idx_q1 = rem

            print(
                f"Mismatch at [{idx_q1}][{j}][{k}][{l}]: C={c_vals[i]}, Rust={rust_vals[i]}"
            )
            diffs += 1
            if diffs > 20:
                print("Too many diffs, stopping")
                break

    if diffs == 0:
        print("Tables match exactly!")
    else:
        print(f"Found {diffs} mismatches.")
