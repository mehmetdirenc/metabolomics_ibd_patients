import os
import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_experiments(folder_path):
    data = {}
    for file in os.listdir(folder_path):
        if file.endswith(".tsv"):
            experiment_name = file.replace(".tsv", "")
            file_path = os.path.join(folder_path, file)
            df = pd.read_csv(file_path, sep="\t")
            expected_columns = ["Compounds", "t.stat", "p.value", "-log10(p)", "FDR", "FoldChange"]
            if not "Compounds" in df.columns:
                # empty_df = pd.DataFrame(columns=expected_columns)
                data[experiment_name] = pd.Series(dtype=float)
            else:
                fold_change_series = np.log2(df.set_index("Compounds")["FoldChange"].replace(0, np.nan))
                data[experiment_name] = fold_change_series
                # data[experiment_name] = np.log2(df.set_index("Compounds")["FoldChange"].replace(0, np.nan))
    return pd.DataFrame(data)

def filter_data(data, exp_type, data_type):
    matching_cols = [col for col in data.columns if col.startswith(exp_type)]
    matching_cols = [col for col in matching_cols if data_type in col]
    if matching_cols:
        return data[matching_cols]  # Drop rows where all selected columns are NaN
    else:
        print(f"No columns found starting with '{exp_type}'.")
        return pd.DataFrame()  # Return an empty DataFrame if no matching columns are found

# def sort_compounds(data):
#     presence_counts = data.notna().sum(axis=1)  # Count non-NaN occurrences across experiments
#     sorted_data = data.loc[presence_counts.sort_values(ascending=False).index]  # Sort by counts
#     return sorted_data

def sort_compounds(data):
    presence_counts = data.notna().sum(axis=1)  # Count non-NaN occurrences across experiments
    mean_fold_change = data.abs().mean(axis=1, skipna=True)  # Calculate mean fold change, ignoring NaNs ABSOLUTE VALUES

    # Create a DataFrame for sorting keys
    sorting_key = pd.DataFrame({
        "presence_counts": presence_counts,
        "mean_fold_change": mean_fold_change
    })

    # Sort by presence counts (descending) and then mean fold change (descending)
    sorting_key = sorting_key.sort_values(by=["presence_counts", "mean_fold_change"], ascending=[False, False], key=abs)

    # Reorder the data based on the sorted key
    sorted_data = data.loc[sorting_key.index]

    return sorted_data

def get_red_blue_extremes_from_sorted(data, top_n=15):
    # Sum of only positive values (reds)
    red_strength = data.clip(lower=0).sum(axis=1, skipna=True)
    # Sum of only negative values (blues)
    blue_strength = data.clip(upper=0).sum(axis=1, skipna=True)  # This will be negative

    # Get top rows by red sum (descending)
    top_red = red_strength.sort_values(ascending=False).head(top_n).index
    # Get top rows by blue sum (ascending, since values are negative)
    top_blue = blue_strength.sort_values().head(top_n).index

    combined_index = list(top_red) + list(top_blue)
    return data.loc[combined_index]

def compare_compounds(df1, df2):
    common_compounds_dict = {}
    all_compounds = []
    # Iterate through the columns in filtered_data
    for col in df1.columns:
        common_col = []
        # Extract the corresponding column name in filtered_data_to_compare
        if "serum" in col:
            compare_col = col.replace('_serum_', '_caecum_')
        elif "caecum" in col:
            compare_col = col.replace('_caecum_', '_serum_')
        else:
            print("problematic data type exit")
            sys.exit()
        # Check if the column exists in filtered_data_to_compare
        if compare_col in df2.columns:
            # Iterate through the rows (compounds) and check for non-NaN values in both columns
            for compound in df1.index:
                if compound in df2[compare_col]:
                    if not pd.isna(df1.loc[compound, col]) and not pd.isna(
                            df2.loc[compound, compare_col]):
                        common_col.append(compound)
                        if compound not in all_compounds:
                            all_compounds.append(compound)
                        # common_col = sorted(list(set(common_col)))
        common_compounds_dict[col] = common_col



    return common_compounds_dict, all_compounds

def create_heatmaps(folder_path, output_folder_path, analysis_type, filtered, data_input_type):
    if filtered:
        output_folder_path = output_folder_path.replace("TO_REPLACE", "filtered") + "_filtered"
        filtered_title = "_filtered"
    else:
        filtered_title = ""
        output_folder_path = output_folder_path.replace("TO_REPLACE", "normal")
  #### FOLD CHANGE DIRECTION IS A/B E.G. CTR_IBD3_CAECUM CITRULLINE = 0.03 A HAS LOWER LEVELS THAN B ####
  #### FOLD CHANGE DIRECTION IS A/B E.G. CTR_IBD3_CAECUM D-GLUCOSAMINE 6-PHOSPHATE = 21.2 A HAS HIGHER LEVELS THAN B ####
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    data_types = ["caecum", "serum"]
    exp_types = ["gno", "muc2_plus", "muc2_minus", "patient1FMT", "patient2FMT", "rescue_pre", "rescue_postplus", "patient_abx"]
    fold_change_data = load_experiments(folder_path)
    fold_change_data = fold_change_data.fillna(np.nan)
    # fold_change_data = fold_change_data.fillna(0)
    # fold_change_data = fold_change_data.loc[~(fold_change_data == 0).all(axis=1)]
    all_dataframes = {}
    for data_type in data_types:
        for exp in exp_types:
            exp_name = exp + "_" + data_type
            if analysis_type == "joined":
                exp_name = exp
                data_type_to_compare = data_types[data_types.index(data_type)-1]
                filtered_data_to_compare = filter_data(fold_change_data, exp, data_type_to_compare)
                # filtered_data_to_compare = filtered_data_to_compare.dropna(how="all")
                filtered_data_to_compare = sort_compounds(filtered_data_to_compare)
                filtered_data_to_compare = filtered_data_to_compare.reindex(sorted(filtered_data_to_compare.columns), axis=1)
                filtered_data = filter_data(fold_change_data, exp, data_type)
                # filtered_data = filtered_data.dropna(how="all")
                filtered_data = sort_compounds(filtered_data)
                common_compounds_dict, all_compounds = compare_compounds(filtered_data, filtered_data_to_compare)
                filtered_data = filtered_data.loc[filtered_data.index.isin(all_compounds)]
                filtered_data = filtered_data.reindex(sorted(filtered_data.columns), axis=1)
                joined_data = filtered_data.join(filtered_data_to_compare, how='left')
                joined_data = joined_data.sort_index(axis=1)
                joined_data = sort_compounds(joined_data)
                joined_data = joined_data.loc[~(joined_data == 0).all(axis=1)]
                joined_data = joined_data.dropna(how="all")
                filtered_data = joined_data
            else:
                heatmap_path = os.path.join(output_folder_path, exp_name)
                filtered_data = filter_data(fold_change_data, exp, data_type)
                filtered_data = sort_compounds(filtered_data)
                filtered_data = filtered_data.reindex(sorted(filtered_data.columns), axis=1)
                filtered_data = filtered_data.loc[~(filtered_data == 0).all(axis=1)]
                filtered_data = filtered_data.dropna(how="all")

            all_dataframes[exp_name] = filtered_data
    # global_min = -10.1
    # global_max = 7
    cases = [["patient_abx_", "patient1FMT_"],
             ["patient2FMT_", "muc2_minus_--hm"],
             ["gno_--ibd1", "muc2_minus_--hm"],
             ["patient_abx_", "muc2_minus_--hm"],
             ["patient1FMT_", "muc2_minus_--hm"]
             ]
    case_frames = {}
    for data_type in data_types:
        for case in cases:
            if analysis_type == "joined":
                data_type = ""
                # case = case[0:-1]
                df_name_1 = case[0].split("--")[0][0:-1] + data_type
                df_name_2 = case[1].split("--")[0][0:-1] + data_type
            else:
                df_name_1 = case[0].split("--")[0] + data_type
                df_name_2 = case[1].split("--")[0] + data_type

            df1 = all_dataframes[df_name_1]
            df2 = all_dataframes[df_name_2]
            joined_df = df1.join(df2, how='inner')
            df_name = df_name_1 + "_" + df_name_2
            all_dataframes[df_name] = joined_df
    # print(all_dataframes.keys())
    global_min = min(all_dataframes[df].min().min() for df in all_dataframes)
    global_max = max(all_dataframes[df].max().max() for df in all_dataframes)
    for exp_name in all_dataframes:
        # if data_input_type == "cmrf" or exp_name.startswith("gno_caecum_m") == False:
        #     continue
        filtered_data = all_dataframes[exp_name]
        # if filtered:
        #     filtered_data = filtered_data.head(40)
        if filtered:
            filtered_data = get_red_blue_extremes_from_sorted(filtered_data, top_n=15)
        heatmap_path = os.path.join(output_folder_path, exp_name)
        title = f"{exp_name}_{data_input_type}_{analysis_type}{filtered_title}"
        if filtered_data.empty:
            fig, ax = plt.subplots(figsize=(20, 15))
            ax.set_xticks([])  # Remove ticks
            ax.set_yticks([])
            ax.set_title(title.upper(), pad=20)
            ax.set_xlabel("Experiments", labelpad=20)
            ax.set_ylabel("Compounds")
            plt.tight_layout()
            plt.savefig(heatmap_path, format='pdf', dpi=600)
            plt.close(fig)
            # plt.show()
            continue
        fig, ax = plt.subplots(figsize=(15, 25))  # Adjust figure size for better layout
        nan_mask = filtered_data.isna()
        sns.heatmap(
            filtered_data,
            cmap="coolwarm",
            linewidths=0.5,
            linecolor="black",
            square=True,  # Ensures equal aspect ratio
            cbar_kws={"label": "Fold Change log2"},
            mask=nan_mask,  # Mask NaN values to keep them white
            ax=ax,
            vmin=global_min,  # Set the minimum value for the color bar
            vmax=global_max,  # Set the maximum value for the color bar
        )
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.set_xlabel("Experiments", labelpad=20)
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        # Label the axes
        plt.title(title.upper(), pad=20)
        plt.xlabel("Experiments")
        plt.ylabel("Compounds")
        plt.xticks(rotation=90)
        # plt.tight_layout()
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        # Save or show the plot
        plt.savefig(heatmap_path, bbox_inches='tight', format='pdf', dpi=600)
        plt.cla()
        plt.close(fig)
        tsv_output_path = heatmap_path + ".tsv"
        filtered_data.to_csv(tsv_output_path, sep="\t", index=True)
        # if exp_name == "gno_caecum":
        #     print("gno_caecum")
        # if "gno_caecum_ours" in title:
        #     plt.show()
        #     print(title)
        # print("x")

if __name__ == '__main__':
    # data_types = ["ours"]
    data_input_types = ["ours", "cmrf" ]
    analysis_types = ["single", "joined"]
    for analysis_type in analysis_types:
        for data_input_type in data_input_types:
            input_folder_path = f"/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics_{data_input_type}"
            output_folder_path = f"/home/direnc/results/mahana/new_metabolomics/TO_REPLACE/metabolomics_heatmaps_{data_input_type}_{analysis_type}"
            create_heatmaps(input_folder_path, output_folder_path, analysis_type, True, data_input_type)
            create_heatmaps(input_folder_path, output_folder_path, analysis_type, False, data_input_type)
    # input_folder_path = "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics_ours"
    # output_folder_path = "/home/direnc/results/mahana/metabolomics_heatmaps_ours_filtered_joined"
    # create_heatmaps(input_folder_path, output_folder_path)