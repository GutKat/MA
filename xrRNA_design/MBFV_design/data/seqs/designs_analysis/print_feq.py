import os

pwd = os.getcwd()

files = [pwd + '/' + i for i in os.listdir(pwd) if i.endswith('.fa')]
files = sorted(files)
for file in files:
	with open(file, 'r') as f:
		print(file)
		for line in f:
			if line.startswith('frequency'):
				print(line)
				break
