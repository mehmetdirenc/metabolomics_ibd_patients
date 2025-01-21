import os











def metaboignite_samplesheet(workdir, samplesheet_path):
    file_list = os.listdir(workdir)
    with open(samplesheet_path, "w") as samplesheet_file:
        header = "sample,level,type,msfile\n"
        samplesheet_file.write(header)
        level = "MS1"
        for file in file_list:
            filepath = os.path.join(workdir, file)
            if file.endswith(".mzML") != True:
                continue
            functional_part = file.split("Neg15_")[1]
            if "blank" in functional_part.lower():
                sample = functional_part.split(".")[0]
                type = "blank"
                samplesheet_file.write("{sample},{level},{type},{msfile}\n".format(sample=sample, level=level, type=type, msfile=filepath))
            else:
                sample = functional_part.split(".")[0]
                type = "disease"
                samplesheet_file.write("{sample},{level},{type},{msfile}\n".format(sample=sample, level=level, type=type, msfile=filepath))





if __name__ == '__main__':
    # workdir = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/data/2024_07_08MS"
    # samplesheet_path = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/data/samplesheet_2024.csv"
    # metaboignite_samplesheet(workdir, samplesheet_path)
    workdir = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/data/2021_CMRF_Mass_Spec_Data"
    samplesheet_path = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/data/samplesheet_2021_CMRF.csv"
    metaboignite_samplesheet(workdir, samplesheet_path)
    pass
