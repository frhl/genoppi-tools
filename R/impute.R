


# replace missing values with randomly sampled values from normal distribution
# with width SD*width and down-shifted MEAN-SD*shift compared to observed sample distribution
# (Perseus, Tyanova et al. 2016)
#def minNormImpute(df,shift,width):
  # df: 'accession_number' column followed by sample columns containing log2-transformed intensity values
  # return same df except with missing values replaced
  
#  temp_df = df.copy()
  
#  temp_df['imputed'] = pd.isnull(temp_df).sum(axis=1) > 0
#  temp_df['imputed'] = temp_df['imputed'].apply(str).str.upper()
  
#  for col in temp_df.columns[1:-1]:
#    mu = np.nanmean(temp_df[col]) - np.nanstd(temp_df[col]) * shift # adjsuted/down-shifted mean (sample mean - SD * shift)
#  sd = np.nanstd(temp_df[col]) * width # adjsuted SD (sample SD * width)
#  temp_df[col] = temp_df[col].apply(lambda x: np.random.normal(mu,sd) if pd.isnull(x) else x)
  
#  return temp_df

impute <- function(df){
  
  
}




