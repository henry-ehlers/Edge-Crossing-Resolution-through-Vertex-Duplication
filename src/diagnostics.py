from pathlib import Path


#
def save_time(times, labels, file_name, directory, overwrite=False):

    # Ensure the length of labels and times are identical
    assert len(times) == len(labels), \
        "Provided times and labels are of different lengths."

    # Create Output Directory if needed
    Path(directory).mkdir(parents=True, exist_ok=True)

    # Open file and append lines to it
    file_path = f"{directory}/{file_name}"
    file = open(file_path, "w") if overwrite else open(file_path, "a")
    for index in range(0, len(times)):
        file.write(f"{labels[index]},{times[index]}\n")
    file.close()
