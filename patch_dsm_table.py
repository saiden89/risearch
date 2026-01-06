dsm_file = "/Users/vzr518/Projects/risearch/src/dsm.rs"
new_table_file = "/Users/vzr518/Projects/risearch/dsm_new_table.rs"

with open(dsm_file, "r") as f:
    dsm_content = f.read()

with open(new_table_file, "r") as f:
    new_table_content = f.read().strip()

# Find the start of DSM_T04_POS
start_marker = "pub const DSM_T04_POS: DsmTable = ["
end_marker = "];"

start_idx = dsm_content.find(start_marker)
if start_idx == -1:
    print("Could not find DSM_T04_POS start")
    exit(1)

# Find the end of the array. Since it is nested, we can't just look for first ];
# We know it ends before DSM_T04_NEG.
next_const = "pub const DSM_T04_NEG: DsmTable = ["
end_idx = dsm_content.find(next_const)

if end_idx == -1:
    # Maybe it's the last one?
    # Let's search for the closing ]; corresponding to the opening [
    # But for now, let's assume valid rust file structure where DSM_T04_NEG follows.
    print("Could not find start of DSM_T04_NEG to delimit end of POS table")
    # Fallback: Find the last ]; before EOF? No.
    # Let's count brackets?
    pass

# Refined extraction:
# We want to replace from start_marker to the ]; before DSM_T04_NEG
# There should be a ]; and some whitespace before DSM_T04_NEG
substring_to_replace_end = dsm_content.rfind("];", 0, end_idx) + 2

if substring_to_replace_end == 1:  # -1 + 2 = 1
    print("Could not find end of table")
    exit(1)

# Construct new content
new_content = (
    dsm_content[:start_idx] + new_table_content + dsm_content[substring_to_replace_end:]
)

with open(dsm_file, "w") as f:
    f.write(new_content)

print("Successfully patched dsm.rs")
