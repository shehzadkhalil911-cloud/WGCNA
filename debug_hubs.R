load('F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/05_hub_genes.RData')
load('F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/01_processed_data_common.RData')

cat('final_hubs length:', length(final_hubs), '\n')
cat('in comm_42568:', sum(final_hubs %in% rownames(comm_42568$expr)), '\n')
cat('in comm_15852:', sum(final_hubs %in% rownames(comm_15852$expr)), '\n')
cat('in comm_29431:', sum(final_hubs %in% rownames(comm_29431$expr)), '\n')
cat('rows 42568:', nrow(comm_42568$expr), '| 15852:', nrow(comm_15852$expr), '| 29431:', nrow(comm_29431$expr), '\n')

# Find which hubs are missing
missing <- final_hubs[!final_hubs %in% rownames(comm_42568$expr)]
cat('missing from 42568:', length(missing), '\n')
print(missing)

# Check expression matrix first few row names
cat('first 10 rownames 42568:\n')
print(head(rownames(comm_42568$expr), 10))
cat('last 10 rownames 42568:\n')
print(tail(rownames(comm_42568$expr), 10))
