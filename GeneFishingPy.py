import numpy as np
import pandas as pd
import utils as utils
from os import listdir
from os.path import isfile, join
from optparse import OptionParser

def main():
	usage = "%prog data_dir_path bait_path output_prefix"
	parser = OptionParser(usage=usage)
	parser.add_option('--alpha', dest='alpha', default=5, 
		help = 'proportion of pool genes to bait genes in each sub iteration: default = %default')
	parser.add_option('--number_of_cluster', dest='number_of_cluster', default=2,
		help = 'number of clusters classified in each sub iteration: default = %default')
	parser.add_option('--way_of_cluster', dest='way_of_cluster', default='Spectral',
		help = 'cluster algorithm, choose between Spectral, K-means and PIC (Power Iteration Clustering)')
	parser.add_option('--affinity_matrix', dest='affinity_matrix', default='Spearman',
		help = 'way ot constructing affinity matrix, choose between Spearman and Cosine')
	parser.add_option('--number_of_iteration', dest='number_of_iteration', default=1000,
		help = 'number of meta iterations: default = %default')
	options, args = parser.parse_args()

	data_dir_path = args[0]
    bait_data_path = args[1]
    output_prefix = args[2]

    alpha = int(options.alpha)
    number_of_cluster = int(options.number_of_cluster)
    way_of_cluster = str(options.way_of_cluster)
    affinity_matrix = str(options.affinity_matrix)
    number_of_iteration = int(options.number_of_iteration)

    if way_of_cluster != 'Spectral' and way_of_cluster != 'PIC' and way_of_cluster != 'K-means':
    	raise ValueError('Wrong clustering method, please choose between Spectral, K-means, and PIC')
    if affinity_matrix != 'Spearman' and affinity_matrix != 'Cosine':
    	raise ValueError('Wrong way of constructing affinity matrix, please choose between Spearman and Cosine')

    #use ABSOLUTE path, you can get it by the command readlink -f file.name
	csv_path = data_dir_path
	csv_files = [f.split('.')[0] for f in listdir(csv_path) if isfile(join(csv_path, f))]
	GTex_Data = {}
	for file in csv_files:
		GTex_Data[file] = pd.read_csv(csv_path + "/"+ file + '.csv', index_col = 0)

	c1_ = pd.read_csv(bait_data_path, header = None).values.tolist()
	c1 = reduce(operator.iconcat, c1_, [])

	for tissues in list(GTex_Data.keys()):
		expr_mat = GTex_Data[tissues]
		bait_genes = list(set(list(expr_mat)).intersection(c1))
		pool_genes = list(set(list(expr_mat)) - set(list(c1)))
		meta_recorder = pd.DataFrame(np.array([np.zeros(len(pool_genes))]), columns = pool_genes)
		for j in range(number_of_iteration):
    		if j % 10 == 0:
        	print('at iteration ' + str(j))
    		utils.gene_fishing(bait_genes, pool_genes, meta_recorder)
			file_name = output_prefix + '/' + tissues + '_CFR.csv'
			meta_recorder.to_csv(file_name)


if __name__ == '__main__':
    main()
