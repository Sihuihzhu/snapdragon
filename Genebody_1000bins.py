import sys
import numpy as np


def print_arr(arr,Chr,endpoint=0,genename= '',order="+"):
	if order == "+":
		arr=list(arr)
		arr.append(endpoint)
		for i in range(len(arr)-1):
			print Chr + "\t" + str(arr[i]) + "\t" + str(arr[i+1]) + "\t" + genename
	else:
		arr=list(arr)
		arr.append(endpoint)
		tmp=[]
		for i in range(len(arr)-1):
			tmp.append('\t'.join([Chr,str(arr[i]),str(arr[i+1]),genename]))
		print "\n".join(list(reversed(tmp)))
	return 0
binsize = 100
bed = file(sys.argv[1])
for line in bed:
	arr=line.strip().split("\t")
	strand = arr[-1]
	genename = arr[3]
	start, end = int(arr[1]), int(arr[2])
	if strand == "+":
		up_interval = np.linspace(start-2000,start,binsize,endpoint=False)
		print_arr(up_interval,arr[0],start,genename,strand)
		gene_body = np.linspace(start,end,binsize,endpoint=False)
		print_arr(gene_body,arr[0],end,genename,strand)
		down_interval = np.linspace(end,end+2000,binsize,endpoint=False)
		print_arr(down_interval,arr[0],end+2000,genename,strand)
		pass
	else:
		up_interval = np.linspace(end,end+2000,binsize,endpoint=False)
		print_arr(up_interval,arr[0],end+2000,genename,strand)
		gene_body = np.linspace(start,end,binsize,endpoint=False)
		print_arr(gene_body,arr[0],end,genename,strand)
		down_interval = np.linspace(start-2000,start,binsize,endpoint=False)
		print_arr(down_interval,arr[0],start,genename,strand)
bed.close()

