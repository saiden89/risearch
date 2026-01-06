import re


def parse_rust_table(content, table_name):
    # Find the const definition
    pattern = r"pub const " + table_name + r": DsmTable = \[(.*?)\];"
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        print(f"Could not find Rust table {table_name}")
        return None

    data_str = match.group(1)
    # Remove comments and whitespace cleanup
    data_str = re.sub(r"//.*", "", data_str)
    # Convert Rust array syntax to Python list syntax
    # Rust uses [ ... ], Python uses [ ... ] - compatible
    # Numbers are just numbers.
    # We need to handle nested brackets carefully or just eval if safe enough
    # Since we trust the file content (it's local source), eval is acceptable for this script
    try:
        # cleanup formatting to make it eval-able
        # Rust might have type annotations or other things, but usually just numbers
        # The DsmTable definition is [[[[i16; 6]; 6]; 6]; 6]
        # In the file it looks like: [ [ [ [-2000, ...], ... ] ... ] ]
        # Only potential issue is trailing commas
        return eval("[" + data_str + "]")
    except Exception as e:
        print(f"Error parsing Rust table {table_name}: {e}")
        return None


def parse_c_table(content, table_name):
    # Find the const definition
    # C syntax: const dsm_t dsm_t04_pos = {\ ... };
    pattern = r"const dsm_t " + table_name + r" = {\\?(.*?)\};"
    match = re.search(pattern, content, re.DOTALL)
    if not match:
        # Try without backslashes in pattern if first failed (some C files might differ)
        pattern = r"const dsm_t " + table_name + r" = {(.*?)};"
        match = re.search(pattern, content, re.DOTALL)
        if not match:
            print(f"Could not find C table {table_name}")
            return None

    data_str = match.group(1)
    # Remove backslashes for line continuation
    data_str = data_str.replace("\\\n", "")
    data_str = data_str.replace("\\", "")
    # Change C braces {} to Python brackets []
    data_str = data_str.replace("{", "[").replace("}", "]")

    try:
        return eval("[" + data_str + "]")
    except Exception as e:
        print(f"Error parsing C table {table_name}: {e}")
        return None


def compare_tables(rust_table, c_table, table_name):
    print(f"Comparing {table_name}")
    print(f"Rust table structure: type={type(rust_table)} len={len(rust_table)}")
    if len(rust_table) > 0:
        print(f"  [0] type={type(rust_table[0])}")
        if isinstance(rust_table[0], list) and len(rust_table[0]) > 0:
            print(f"  [0][0] type={type(rust_table[0][0])}")

    print(f"C table structure: type={type(c_table)} len={len(c_table)}")

    mismatches = 0

    # 6x6x6x6
    bases = ["GAP", "A", "G", "C", "U", "N"]

    for q1 in range(6):
        for q2 in range(6):
            for t1 in range(6):
                for t2 in range(6):
                    r_val = rust_table[q1][q2][t1][t2]
                    try:
                        c_val = c_table[q1][q2][t1][t2]
                    except IndexError:
                        print(f"IndexError in C table at [{q1}][{q2}][{t1}][{t2}]")
                        return

                    if r_val != c_val:
                        mismatches += 1
                        print(
                            f"MISMATCH in {table_name} at Q:{bases[q1]}{bases[q2]} T:{bases[t1]}{bases[t2]} ([{q1}][{q2}][{t1}][{t2}]) :: Rust={r_val}, C={c_val}"
                        )

    if mismatches == 0:
        print(f"SUCCESS: {table_name} matches perfectly.")
    else:
        print(f"FAILURE: {table_name} had {mismatches} mismatches.")


def main():
    rust_path = "src/dsm.rs"
    c_t04_path = "legacy_c/RIsearch2/src/dsm-t04-guugle.c"
    c_ext_path = "legacy_c/RIsearch2/src/dsm-extend-guugle.c"

    try:
        with open(rust_path, "r") as f:
            rust_content = f.read()
        with open(c_t04_path, "r") as f:
            c_t04_content = f.read()
        with open(c_ext_path, "r") as f:
            c_ext_content = f.read()
    except FileNotFoundError as e:
        print(f"Error opening file: {e}")
        return

    # Mappings
    # Rust Name -> (C Content, C Name)
    tasks = [
        ("DSM_T04_POS", c_t04_content, "dsm_t04_pos"),
        ("DSM_T04_NEG", c_t04_content, "dsm_t04_neg"),
        ("DSM_EXTEND_POS", c_ext_content, "dsm_extend_pos"),
        ("DSM_EXTEND_NEG", c_ext_content, "dsm_extend_neg"),
    ]

    for r_name, c_content, c_name in tasks:
        r_table = parse_rust_table(rust_content, r_name)
        c_table = parse_c_table(c_content, c_name)

        if r_table and c_table:
            # r_table and c_table are already the [6][6][6][6] structure
            compare_tables(r_table, c_table, r_name)


if __name__ == "__main__":
    main()
