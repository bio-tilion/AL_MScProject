import pandas as pd


class ExpressionMatrix:
    def __init__(self, df: pd.DataFrame, metadata: pd.DataFrame):
        if df.iloc[:,1:].shape[1] == metadata.shape[0]:
            self.expression = df
            self.metadata = metadata
        else:
            raise Exception("Metadata dataframe shape is not congruent with expression dataframe shape.")
    
    def filter_by_metadata(self, metadata_df):
        # get selection of column labels (first column) from metadata
        selection = metadata_df.iloc[:, 0]

        # filter expression dataset
        df_out = self.expression.loc[:, selection]

        # add back (first) column with gene identifiers
        gene_column = self.expression.iloc[:, 0]
        df_out.insert(0, gene_column.name, gene_column)

        return ExpressionMatrix(df_out, metadata_df)
    
    def filter_by_metadata_value(self, field: str, value: str):
        # get metadata table
        metadata_out = self.metadata[self.metadata[field] == value]
        
        return self.filter_by_metadata(metadata_df=metadata_out)
    
    def filter_by_gene(self, subset: list):
        # filter for subset list (assuming gene identifiers are in the first column)
        df_out = self.expression[self.expression.iloc[:, 0].isin(subset)]

        return ExpressionMatrix(df_out, self.metadata)
    
    def save_expression(self, file_path: str, metadata: bool = False):
        if metadata:
            path_metadata = file_path + "_metadata.csv"
        path_expression = file_path + ".txt"
        
        try:
            txt_metadata = ""
            self.expression.to_csv(path_expression, index=False, sep="\t")
            txt_expression = f"Expression matrix saved to {path_expression}"
            if metadata:
                self.metadata.to_csv(path_metadata, index=False)
                txt_metadata = f"\nMetadata matrix saved to {path_metadata}"
            txt_out = f"{txt_expression}{txt_metadata}"
        except:
            txt_out = "Something went wrong"
        
        print(txt_out)
        return None
    
    def has_unique_time_measurments(self, time_column: str = "time") -> bool:
        return len(self.metadata[time_column]) == len(self.metadata[time_column].unique())

    def convert_labels_to_time(self, label_column: str = "label", time_column: str = "time") -> None:
        if self.has_unique_time_measurments(time_column = time_column):
            column_labels = dict()
            for label, timepoint in zip(self.metadata[label_column], self.metadata[time_column]):
                column_labels[label] = timepoint
            
            # copy first to get rid of warning
            metadata_temp = self.metadata.copy()
            metadata_temp[label_column] = metadata_temp[time_column]

            self.metadata = metadata_temp
            self.expression.rename(columns=column_labels, inplace=True)
            return None
        else:
            raise ValueError("The dataset doesn't have a single measurement per timepoint")
