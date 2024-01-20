mkdir -p Save/fits
mkdir -p Save/pdep_post
mkdir -p Data/analyzed
mkdir -p Data/curia_data
mkdir -p Data/track2_20220404

aws s3 cp s3://akassler-bcbs-curia-causal/2023_12_10_df_merged_raw_practice_level_10x.zip .
unzip 2023_12_10_df_merged_raw_practice_level_10x.zip
rm -f 2023_12_10_df_merged_raw_practice_level_10x.zip
mv df_merged_raw_practice_level ./Data/curia_data/df_merged_raw_practice_level_10x
aws s3 cp s3://akassler-bcbs-curia-causal/dataset_nums.csv ./Data/curia_data/
aws s3 cp s3://akassler-bcbs-curia-causal/track2_20220404.zip ./Data/
unzip ./Data/track2_20220404.zip
rm -f ./Data/track2_20220404.zip
mv practice Data/track2_20220404
mv practice_year Data/track2_20220404
aws s3 cp s3://akassler-bcbs-curia-causal/sdY.rds .
aws s3 cp s3://akassler-bcbs-curia-causal/run_all_loops.R .

yum install -y tmux

#aws s3 cp s3://akassler-bcbs-curia-causal/fit_bart_functional_cloud.R .
#sudo yum update -y
#sudo yum install -y R
