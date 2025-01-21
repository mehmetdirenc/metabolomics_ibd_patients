import glob
import os, io
import sys
from pyclbr import Class

import pandas as pd

def write_sheets(excel_doc, output_folder):
    for sheet_name, df in pd.read_excel(excel_doc, sheet_name=None).items():
        df.to_csv(os.path.join(output_folder, f'{sheet_name}.csv'), index=False)


def import_mass_finder_results(mass_finder_mztab):
    filtered_lines = []
    with open(mass_finder_mztab, "r") as file:
        for line in file:
            # Skip lines starting with "MTD" or empty lines
            if not line.startswith("MTD") and line.strip():
                filtered_lines.append(line)
    mass_df = pd.read_csv(io.StringIO("".join(filtered_lines)), sep="\t")
    return mass_df


def select_row_for_description(description, filtered_df, cp_df):
    # if description == "uracil":
    #     print("a")
    # Check if the description exists in cp_df
    cp_row = cp_df[cp_df["compound"].str.lower() == description.lower()]
    # y = filtered_df["description"]
    # x = description in filtered_df["description"].str.lower()
    if not cp_row.empty:
        # Get medRt_seconds for the compound in cp_df
        med_rt = cp_row.iloc[0]["medRt_seconds"]
        matching_rows = filtered_df[filtered_df["description"].str.lower() == description.lower()]

        if not matching_rows.empty:
            # Find the row with the closest retention_time to med_rt
            closest_idx = (matching_rows["retention_time"] - med_rt).abs().idxmin()
            closest_row = matching_rows.loc[closest_idx]
            return closest_row
        # Find the closest retention_time in filtered_df
        # closest_row = filtered_df.loc[(filtered_df["description"].str.lower() == description.lower()) &
        #                               (filtered_df["retention_time"] - med_rt).abs().idxmin()]
            return closest_row
    else:
        # If not in cp_df, select the row with the highest abundance
        highest_abundance_row = filtered_df.loc[filtered_df["description"].str.lower() == description.lower()].nlargest(
            1, "smallmolecule_abundance_study_variable[1]"
        )
        return highest_abundance_row.iloc[0] if not highest_abundance_row.empty else None


def filter_group_compounds(df, cp_df):
    filtered_df = df[(df["retention_time"] >= 0.5) & (df["calc_mass_to_charge"] >= 5) & (df["smallmolecule_abundance_study_variable[1]"] >= 10000)]
    unique_descriptions = filtered_df["description"].str.lower().unique()

    result_rows = [
        select_row_for_description(description, filtered_df, cp_df)
        for description in unique_descriptions
    ]

    # Create the final DataFrame with the selected rows
    result_df = pd.DataFrame([row for row in result_rows if row is not None])
    # grouped_df = filtered_df.groupby("description", as_index=False)["smallmolecule_abundance_study_variable[1]"].sum()
    return result_df

def mass_finder_files_dict(year_mass_finder, identifier, openms_results_folder):
    if year_mass_finder == "2024":
        mass_finder_results_path = os.path.join(openms_results_folder, "2024_mass_finder")
    elif year_mass_finder == "2021":
        mass_finder_results_path = os.path.join(openms_results_folder, "2021_mass_finder")
    else:
        print("problematic year %s"%year_mass_finder)
        sys.exit()
    mztab_file_dict = {}
    mztab_files = glob.glob(os.path.join(mass_finder_results_path, "*.mzTab"))
    for mztab_file in mztab_files:
        if year_mass_finder == "2021":
            if "Sample" not in mztab_file:
                continue
            sample_no = mztab_file.split("_Sample")[1].split(".")[0]
            mztab_file_dict[sample_no] = mztab_file
        elif identifier in mztab_file:
            sample_no = mztab_file.split(identifier)[1].split("_")[0]
            if not sample_no.isdigit():
                print(Class(sample_no))
                print("mztab_file sample no issue %s"%mztab_file)
                sys.exit(1)
            else:
                mztab_file_dict[sample_no] = mztab_file
    return mztab_file_dict


def compounds_expected_rt(rt_compound_filepath):
    cp_df = pd.read_csv(rt_compound_filepath, sep=",")
    return cp_df

def get_sample_filename(oif_df, identifier, openms_results_folder, sample_no, year):
    mztab_file_dict = mass_finder_files_dict(year, identifier, openms_results_folder)
    mztab_file = mztab_file_dict[sample_no]
    return mztab_file



# def analyze_mass_finder(mass_finder_mztab, old_intensities_folder, new_intensities_folder, openms_results_folder):
#     mass_df = import_mass_finder_results(mass_finder_mztab)
#     new_df = filter_group_compounds(mass_df)
#     # adapt_intensities(old_intensities_folder, new_intensities_folder, new_df, openms_results_folder_2021, openms_results_folder_2024)
#     return


def import_cp_rt(rt_compound_filepath_2024, rt_compound_filepath_2021):
    cp_df_2021 = pd.read_csv(rt_compound_filepath_2021, sep=",")
    cp_df_2024 = pd.read_csv(rt_compound_filepath_2024, sep=",")
    cp_df_2021["medRt_seconds"] = cp_df_2021["medRt"] * 60
    cp_df_2024["medRt_seconds"] = cp_df_2024["medRt"] * 60
    return cp_df_2021, cp_df_2024

def adapt_intensities(old_intensities_folder, new_intensities_folder, openms_results_folder, rt_compound_filepath_2024, rt_compound_filepath_2021):
    cp_df_2021, cp_df_2024 = import_cp_rt(rt_compound_filepath_2024, rt_compound_filepath_2021)
    analyzed_files = []
    if not os.path.exists(new_intensities_folder):
        os.mkdir(new_intensities_folder)
    old_intensity_files = glob.glob(os.path.join(old_intensities_folder, "*.csv"))
    for oif in old_intensity_files:
        basename_oif = os.path.basename(oif)
        res_csv_path = os.path.join(new_intensities_folder, basename_oif)
        print(basename_oif)
        if os.path.exists(res_csv_path):
            continue
        # Initialize a dictionary to store data
        data_dict = {}

        # Collect unique compounds and populate data_dict
        all_compounds = set()
        if "serum" in oif:
            identifier = "_S_"
        elif "caecum" in oif:
            identifier = "_C_"
        else:
            continue
        oif_df = pd.read_csv(oif, sep=",", header=[0])
        file_cols = oif_df.columns.tolist()
        labels = oif_df.loc[0, :].values.tolist()
        year = ""
        for sample in file_cols:
            if sample == "Sample":
                continue
            if "sample" in sample.lower():
                year = "2021"
            elif "s_" in sample.lower() or "c_" in sample.lower():
                year = "2024"
            else:
                print("problematic value%s" % sample)
                sys.exit()
            if year == "2021":
                sample_no = sample.strip("Sample")
            else:
                if identifier == "_S_":
                    sample_no = sample.split("S_")[1].split("_")[0]
                else:
                    sample_no = sample.split("C_")[1].split("_")[0]
            # except:
            #     print("issue")
            mztab_file = get_sample_filename(oif_df, identifier, openms_results_folder, sample_no, year)
            mass_df = import_mass_finder_results(mztab_file)
            # mass_df.drop(mass_df.tail(5).index,
            #         inplace=True)
            if year == "2024":
                new_df = filter_group_compounds(mass_df, cp_df_2024)
            else:
                new_df = filter_group_compounds(mass_df, cp_df_2021)
            compounds = new_df['description']
            abundances = new_df['smallmolecule_abundance_study_variable[1]']
            # Map compounds to their abundances for this sample
            compound_dict = dict(zip(compounds, abundances))
            data_dict[sample] = compound_dict
            all_compounds.update(compounds)
            # print("asd")
        # pass
        final_df = pd.DataFrame(index=list(all_compounds))

        # Add abundance data for each sample
        for sample in file_cols:
            if sample == "Sample":
                continue
            # Map abundances to the final dataframe, filling missing values with 0
            final_df[sample] = final_df.index.map(data_dict.get(sample, {}))
            final_df.dropna(subset=[sample], inplace=True)

        # Add the labels row
        # final_df.loc['Label'] = labels
        # final_df.loc['Sample'] = file_cols
        final_df.reset_index(inplace=True)
        final_df.rename(columns={'index': 'Compound'}, inplace=True)
        # Create a new DataFrame for the 'Label' and 'Sample' rows
        sample_row = pd.DataFrame([file_cols], columns=final_df.columns)
        labels_row = pd.DataFrame([labels], columns=final_df.columns)

        # Append 'Label' and 'Sample' rows to the beginning of the final dataframe
        final_df = pd.concat([labels_row, sample_row, final_df])

        # Remove quotes from the first column (if any)
        # Reset index to move 'Compound' to the first column
        final_df['Compound'] = final_df['Compound'].str.replace(',', '_', regex=False)
        final_df['Compound'] = final_df['Compound'].str.replace('"', '', regex=False)
        # Save to CSV
        final_df.to_csv(res_csv_path, index=False, header=False)
    return

# def adjust_pos


if __name__ == '__main__':
    # mztab_file = "/home/direnc/inputs/mahana/ms_ms/2024_07_08_MS_QEB-HILICNeg15_C_11_cal8pre.mzTab"
    # mztab_file = "/mnt/lustre/home/mager/magmu818/results/mahana/ms/2024_results/compounds/2024_07_08_MS_QEB-HILICNeg15_C_11_cal8pre.mzTab"
    # excel_doc = r"/home/direnc/inputs/mahana/ms_ms/metabolomic_explanations_restructured.xlsx"
    # output_folder = '/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/'
    # write_sheets(excel_doc, output_folder)
    ### lOCAL ###
    # rt_compound_filepath_2024 = "/home/direnc/inputs/mahana/ms_ms/2024_compounds.csv"
    # rt_compound_filepath_2021 = "/home/direnc/inputs/mahana/ms_ms/2021_compounds.csv"
    # openms_results_folder_2021 = "/home/direnc/inputs/mahana/ms_ms/2021_mass_finder"
    # openms_results_folder_2024 = "/home/direnc/inputs/mahana/ms_ms/2024_mass_finder"
    # openms_results_folder = "/home/direnc/inputs/mahana/ms_ms/"
    # old_intensities_folder = "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics_cmrf/"
    # new_intensities_folder = "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics_ours/"
    ### SERVER ###
    rt_compound_filepath_2024 = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/2024_compounds.csv"
    rt_compound_filepath_2021 = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/2021_compounds.csv"
    openms_results_folder_2021 = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/2021_mass_finder"
    openms_results_folder_2024 = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/2024_mass_finder"
    openms_results_folder = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/"
    old_intensities_folder = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/metabolomics_statistics_cmrf/"
    new_intensities_folder = "/mnt/lustre/home/mager/magmu818/inputs/mahana/ms/ms_ms_analysis_files/metabolomics_statistics_ours/"
    adapt_intensities(old_intensities_folder, new_intensities_folder, openms_results_folder, rt_compound_filepath_2024, rt_compound_filepath_2021)
    # analyze_mass_finder(mztab_file, old_intensities_folder, new_intensities_folder, openms_results_folder)
    pass