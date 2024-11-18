#import libraries + data
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st

#######################configure streamlit#########################
def main():
    st.title('Cas9 Comparison')

    #file upload
    st.sidebar.header('Upload Data File')
    uploaded_file = st.sidebar.file_uploader('Upload an Excel file', type=['xlsx'])

    if uploaded_file is None:
        st.warning("Please upload an Excel file to begin.")
        return

    if uploaded_file is not None:
        st.subheader('Uploaded Data')
        df = import_data_clean(uploaded_file)
    
    #prepare data
    try:
        df, variant_cols = prepare_dataframe(df)
    except ValueError as e:
        st.error(f'Error: {e}')
        return
    
    #create subsets
    df_on, df_off = create_subsets(df, variant_cols)

    #show subset data
    st.subheader('On-Target-Data')
    st.write(df_on.head())
    st.subheader('Off-Target-Data')
    st.write(df_off.head())

    #fidelity + efficiency calculations
    st.subheader('Fidelity and Efficiency Results')
    fidelity_efficiency_variants = fid_eff_variants(df_on, df_off, variant_cols)
    fidelity_efficiency_gRNAs = fid_eff_gRNAs(df_on, df_off, variant_cols)

    st.write('Fidelity and Efficiency for variants:')
    st.write(fidelity_efficiency_variants)

    st.write('Fidelity and Efficiency for gRNAs:')
    st.write(fidelity_efficiency_gRNAs)

    #plot results
    st.subheader('Fidelity vs Efficiency Plot for Variants')
    fig = plot_fid_eff_variants(fidelity_efficiency_variants)
    st.pyplot(fig)

    #option to download results
    st.subheader('Download Results')
    st.download_button(
        label="Download On-Target Data",
        data = df_on.to_csv(index=False),
        file_name='on_target_data.csv',
        mime='text/csv'
    )
    st.download_button(
        label="Download Off-Target Data",
        data=df_off.to_csv(index=False),
        file_name="off_target_data.csv",
        mime="text/csv"
    )
    st.download_button(
        label="Download Variant Fidelity and Efficiency",
        data=fidelity_efficiency_variants.to_csv(index=False),
        file_name="fidelity_efficiency_variants.csv",
        mime="text/csv"
        )
    st.download_button(
        label="Download gRNA Fidelity and Efficiency",
        data=fidelity_efficiency_gRNAs.to_csv(index=False),
        file_name="fidelity_efficiency_gRNAs.csv",
        mime="text/csv"
        )

    



############################utility functions####################
#import data
def import_data_clean(filename):

    #Dynamically detect the start of the dataframe in an Excel file.
    # Read the first few rows to detect the header
    preview = pd.read_excel(filename, nrows=10, header=None)
    # Find the row where the table starts
    header_row_index = None
    for index, row in preview.iterrows():
        if 'MMs' in row.values:
            header_row_index = index
            break
    
    if header_row_index is None:
        raise ValueError("Could not locate the start of the table dynamically.")

    df = pd.read_excel(filename, skiprows=header_row_index,header=0)
    return df

#identify relevant columns (gRNA, MMs, variants ('Cas9'))
def prepare_dataframe(df):
    df.columns = df.columns.astype(str)
    gRNA_col = [col for col in df.columns if 'grna' in col.lower()]
    
    # Check if a gRNA column is found
    if len(gRNA_col) == 0:
        raise ValueError("No gRNA column found in the DataFrame.")
    elif len(gRNA_col) > 1:
        raise ValueError("Multiple gRNA columns found: {}".format(gRNA_col))
    
    gRNA_col = gRNA_col[0]  # Use the first found gRNA column
    df.rename(columns={gRNA_col: 'gRNA'}, inplace=True)

    MMs_col = 'MMs'
    variant_cols = [col for col in df.columns if col not in [gRNA_col, MMs_col] and 'cas' in col.lower()]

    #numeric data for variant columns
    for col in variant_cols + [MMs_col]:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

    #total indels column by summing across all set variant columns
    df['total indels'] = df[variant_cols].sum(axis=1)
    return df, variant_cols




# Subsets for on-target and off-target data
def create_subsets(dataframe, variant_cols):
    # On-target data where MMs is 0
    df_on = dataframe[dataframe["MMs"] == 0].reset_index(drop=True)

    # Off-target data, copy the original DataFrame
    df_off = dataframe.copy()

    #numeric columns subset
    numeric_cols = df_off.select_dtypes(include=[np.number]).columns.tolist()

    # Set variant columns to 0 where MMs is 0
    for col in variant_cols:
        df_off[col] = np.where(df_off["MMs"] == 0, 0, df_off[col])
    
    #set total indels to 0
    for col in numeric_cols:
        df_off[col] = np.where(df_off['MMs'] == 0, 0, df_off[col])

    

    return df_on, df_off




#fidelity and efficiency for variants
def fid_eff_variants(ondata, offdata, variant_cols):
    sum_on_variants = ondata[variant_cols].sum()
    sum_off_variants = offdata[variant_cols].sum()
    total_on_count = ondata['total indels'].sum()
    average_on_count = total_on_count / len(variant_cols)


    fid_variants = sum_on_variants / (sum_on_variants + sum_off_variants)
    eff_variants = sum_on_variants / average_on_count

    #results df (columns(summed on and off-target indels + average on-target count): to retrace fidelity and efficiency calculations)
    fid_eff_variants_df = pd.DataFrame({
        'variant': variant_cols,
        'summed on-target indels': sum_on_variants.values,
        'summed off-target indels': sum_off_variants.values,
        'fidelity variant': fid_variants.values,
        'efficiency variant': eff_variants.values,
        'average on-target count per variant': average_on_count
    })
    return fid_eff_variants_df



#fidelity and efficiency for gRNAs
def fid_eff_gRNAs(ondata, offdata, variant_cols):
    sum_on_gRNAs = ondata['total indels']  #total indel count for each gRNA
    summed_groups = offdata.groupby('gRNA', sort=False)[variant_cols].agg('sum').reset_index()
    sum_off_gRNAs = summed_groups.iloc[:, 1:].sum(axis=1)
    gRNA_names = summed_groups['gRNA']
    
    #average on-target indel count per gRNA
    average_on_count_gRNA = sum_on_gRNAs.sum() / len(sum_on_gRNAs)

    # Calculate fidelity and efficiency for gRNAs
    fid_gRNAs = sum_on_gRNAs / (sum_on_gRNAs + sum_off_gRNAs)
    eff_gRNAs = sum_on_gRNAs / average_on_count_gRNA

    fid_eff_gRNAs_df = pd.DataFrame({
        'gRNA': gRNA_names,
        'summed on-target indels': sum_on_gRNAs.values,
        'summed off-target indels': sum_off_gRNAs.values,
        'fidelity gRNA': fid_gRNAs.values,
        'efficiency gRNA': eff_gRNAs.values,
        'average on-target count per gRNA': average_on_count_gRNA
    })
    return fid_eff_gRNAs_df



#plot fidelity and efficiency for variants
def plot_fid_eff_variants(fid_eff_variants_df):
    #convert fidelity and efficiency to percentages
    fid_eff_variants_df['fidelity variant (%)'] = fid_eff_variants_df['fidelity variant'] * 100
    fid_eff_variants_df['efficiency variant (%)'] = fid_eff_variants_df['efficiency variant'] * 100

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.scatterplot(
        x='fidelity variant (%)', 
        y='efficiency variant (%)', 
        data=fid_eff_variants_df, 
        hue='variant', 
        palette='tab10', 
        s=100, 
        legend=False,
        ax=ax
    )

    #labels for each variant
    for i, row in fid_eff_variants_df.iterrows():
        ax.text(
            row['fidelity variant (%)'] + 1, 
            row['efficiency variant (%)'] - 6, 
            row['variant'], 
            horizontalalignment='right', 
            size='medium', 
            color='black', 
            weight='semibold'
        )

    #dynamically adjust the limits for the x and y axes based on the data
    x_min = fid_eff_variants_df['fidelity variant (%)'].min() - 5
    x_max = fid_eff_variants_df['fidelity variant (%)'].max() + 5
    y_min = fid_eff_variants_df['efficiency variant (%)'].min() - 10
    y_max = fid_eff_variants_df['efficiency variant (%)'].max() + 10

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    #dynamically set xticks and yticks (np.arange takes the )
    x_ticks = np.arange(np.floor(x_min/10)*10, np.ceil(x_max/10)*10, 10)
    y_ticks = np.arange(np.floor(y_min/10)*10, np.ceil(y_max/10)*10, 10)

    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)

    #labels and title
    ax.set_xlabel("Fidelity (%)")
    ax.set_ylabel("Efficiency (%)")
    ax.set_title("Fidelity and Efficiency Scores for Variants")
    ax.grid(True, linestyle='--', alpha=0.7)

    return fig

if __name__ == "__main__":
    main()





