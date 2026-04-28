import os
import re


def read_alpha_water(filepath):
    """Reads alpha.water and returns the list of scalar values."""
    with open(filepath, "r") as f:
        content = f.read()

    # Find the nonuniform List<scalar> block
    # Pattern: nonuniform List<scalar> \n <size> \n ( ... )
    match = re.search(
        r"nonuniform List<scalar>\s*\n\s*(\d+)\s*\n\s*\(([\s\S]*?)\)", content
    )
    if not match:
        print("Error: Could not find nonuniform List<scalar> in alpha.water")
        return None, None

    size = int(match.group(1))
    data_str = match.group(2)

    # Parse values
    values = []
    for val in data_str.split():
        values.append(float(val))

    if len(values) != size:
        print(f"Error: Expected {size} values, found {len(values)}")
        return None, None

    return size, values


def create_co2_field(alpha_values):
    """Creates C_O2 values based on alpha.water."""
    co2_values = []
    for alpha in alpha_values:
        if alpha >= 0.5:
            co2_values.append(0.0)  # Water
        else:
            co2_values.append(9.0)  # Air
    return co2_values


def write_co2_file(template_path, output_path, size, co2_values):
    """Reads template C_O2 and writes new file with updated internalField."""
    with open(template_path, "r") as f:
        content = f.read()

    # Construct the new internalField string
    new_data_str = "\n".join([str(val) for val in co2_values])
    new_internal_field = (
        f"internalField   nonuniform List<scalar> \n{size}\n(\n{new_data_str}\n);"
    )

    # Replace the existing internalField
    # It might be 'uniform ...' or 'nonuniform ...'
    # strict regex to replace specifically the internalField entry

    # Pattern to match: internalField <anything until semicolon>;
    # Using non-greedy match including newlines
    new_content = re.sub(
        r"internalField\s+[\s\S]*?;", new_internal_field, content, count=1
    )

    with open(output_path, "w") as f:
        f.write(new_content)


def main():
    alpha_path = "0/alpha.water"
    co2_template = "0/C_O2"  # We will overwrite this or use it as template

    print("Reading alpha.water...")
    size, alpha_values = read_alpha_water(alpha_path)
    if not size:
        return

    print("Generating C_O2 field...")
    co2_values = create_co2_field(alpha_values)

    print(f"Writing updated {co2_template}...")
    write_co2_file(co2_template, co2_template, size, co2_values)
    print("Done.")


if __name__ == "__main__":
    main()
