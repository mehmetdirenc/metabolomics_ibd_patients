import os
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
    mean_fold_change = data.mean(axis=1, skipna=True)  # Calculate mean fold change, ignoring NaNs

    # Create a DataFrame for sorting keys
    sorting_key = pd.DataFrame({
        "presence_counts": presence_counts,
        "mean_fold_change": mean_fold_change
    })

    # Sort by presence counts (descending) and then mean fold change (descending)
    sorting_key = sorting_key.sort_values(by=["presence_counts", "mean_fold_change"], ascending=[False, False])

    # Reorder the data based on the sorted key
    sorted_data = data.loc[sorting_key.index]

    return sorted_data


def create_heatmaps(folder_path, output_folder_path):
    data_types = ["caecum", "serum"]
    exp_types = ["gno", "muc2_plus", "muc2_minus", "patient1FMT", "patient2FMT", "rescue_pre", "rescue_postplus", "patient_abx"]
    fold_change_data = load_experiments(folder_path)
    fold_change_data = fold_change_data.fillna(np.nan)
    for data_type in data_types:
        for exp in exp_types:
            exp_name = exp + "_" + data_type
            heatmap_path = os.path.join(output_folder_path, exp_name)
            filtered_data = filter_data(fold_change_data, exp, data_type)
            filtered_data = filtered_data.dropna(how="all")
            filtered_data = sort_compounds(filtered_data)
            filtered_data = filtered_data.reindex(sorted(filtered_data.columns), axis=1)
            # print(f"Filtered data for experiment type: {exp}")
            # print(filtered_data)
        # for exp_type in exp_types:
        #     fold_change_data = load_experiments(folder_path, exp_type)
        #     fold_change_data = fold_change_data.fillna(np.nan)  # Keep NaN for missing values

        # Step 3: Plot the heatmap
        #     plt.figure(figsize=(40, 40 ))  # Adjust size to ensure square cells
            if filtered_data.empty:
                fig, ax = plt.subplots(figsize=(20, 15))
                ax.set_xticks([])  # Remove ticks
                ax.set_yticks([])
                ax.set_title(exp_name.upper(), pad=20)
                ax.set_xlabel("Experiments", labelpad=20)
                ax.set_ylabel("Compounds")
                plt.tight_layout()
                plt.savefig(heatmap_path)
                # plt.show()
                continue
            fig, ax = plt.subplots(figsize=(7, 35))  # Adjust figure size for better layout
            nan_mask = filtered_data.isna()
            sns.heatmap(
                filtered_data,
                cmap="coolwarm",
                linewidths=0.5,
                linecolor="black",
                square=True,  # Ensures equal aspect ratio
                cbar_kws={"label": "Fold Change"},
                mask=nan_mask,  # Mask NaN values to keep them white
                ax=ax,
            )
            for i in range(filtered_data.shape[0]):  # Iterate over rows
                for j in range(filtered_data.shape[1]):  # Iterate over columns
                    if nan_mask.iloc[i, j]:  # If the value is NaN
                        ax.text(
                            j + 0.5, i + 0.5, "X",  # Position the text in the center of the box
                            ha="center", va="center", color="black", fontsize=10, fontweight="bold"
                        )
            ax.xaxis.set_label_position("top")
            ax.xaxis.tick_top()
            ax.set_xlabel("Experiments", labelpad=20)  # Add padding for clarity
            # Remove any unnecessary extra bars
            # plt.tight_layout()  # Ensure proper layout without extra gaps

            # Add hatch lines for missing values
            # missing_mask = fold_change_data.isna()
            # sns.heatmap(
            #     missing_mask,
            #     cmap=sns.color_palette(["lightgrey"], as_cmap=True),
            #     annot=False,
            #     linewidths=0.5,
            #     linecolor='black',
            #     alpha=0.2,  # Adjust transparency
            # )

            ax.spines['top'].set_visible(True)
            ax.spines['right'].set_visible(True)
            ax.spines['bottom'].set_visible(True)
            ax.spines['left'].set_visible(True)
            # Label the axes
            plt.title(exp_name.upper(), pad=20)
            plt.xlabel("Experiments")
            plt.ylabel("Compounds")
            plt.xticks(rotation=90)
            plt.tight_layout()
            # Save or show the plot
            # plt.savefig(heatmap_path)
            plt.show()
            print("x")

if __name__ == '__main__':
    input_folder_path = "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics"
    output_folder_path = "/home/direnc/results/mahana/metabolomics_heatmaps"
    create_heatmaps(input_folder_path, output_folder_path)