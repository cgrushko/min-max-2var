#!/usr/bin/python
import pprint, os, json
pp = pprint.PrettyPrinter(indent=4)

import subprocess, re, texttable, sys, numpy
from collections import defaultdict

if len(sys.argv) < 2:
	print '(OPTIONAL) export TEXTTABLE_WIDTH=150'
	print './Usage: aggregate_benchmark.results.py FILE1 FILE2 ...'
	quit()

results = defaultdict(lambda : defaultdict(list))
results_table = list()
test_names = set([])
solver_names = set(['hough'])

re1 = re.compile(r'Hough total time \(ms\): ([\d.]+)')  # Hough total time (ms): 81932.702000
re2 = re.compile(r'([\w \(\),]+?) total time \(ms\): ([\d.]+)')  # lpsolve55 (etaPFI) total time (ms): 81932.702000

sys.argv.pop(0)
for arg in sys.argv:
	times_text = subprocess.check_output(['sh', './extract_times_from_benchmarks.sh', arg])
	times_text = times_text.split('\n')
	for i in range(0, len(times_text)-1, 3):
		test_name = times_text[i]
		hough_time = float(re1.match(times_text[i+1]).groups()[0]) / 1000
		external_solver_name, external_solver_time = re2.match(times_text[i+2]).groups()
		external_solver_time = float(external_solver_time) / 1000
		results[test_name]['hough'].append(hough_time)
		results[test_name][external_solver_name].append(external_solver_time)
		solver_names.add(external_solver_name)
		test_names.add(test_name)

#############

if 'TEXTTABLE_WIDTH' in os.environ:
	texttable_width = int(os.environ['TEXTTABLE_WIDTH'])
else:
	texttable_width = 140
table = texttable.Texttable(texttable_width)
row = list(test_names)
row.insert(0, '')
table.add_row(row)
hough_avg_time = dict()
for test_name in test_names:
	hough_avg_time[test_name] = numpy.mean(results[test_name]['hough'])
for solver_name in solver_names:
	row = [solver_name]
	for test_name in test_names:
		row.append(', '.join(map('{0:.2f}'.format, results[test_name][solver_name])))
	table.add_row(row)		

print '(In seconds)'
print table.draw() + "\n"

#############

table = texttable.Texttable(140)
row = list(test_names)
row.insert(0, '')
table.add_row(row)
for solver_name in solver_names:
	if solver_name == 'hough':
		continue
	row = [solver_name]
	for test_name in test_names:
		external_solver_avg_time = numpy.mean(results[test_name][solver_name])
		row.append('x {0:.2f}'.format(external_solver_avg_time / hough_avg_time[test_name]))
	table.add_row(row)		

print '(Times slower)'
print table.draw() + "\n"

############# Output JS for highcharts (see chart_template.js)

def round_to_2_digits(x):
	return round(x * 100) / 100

test_names = sorted(test_names)
# print str(test_names)
result_list = []
for solver_name in solver_names:
	row = {'name': solver_name, 'data':[]}
	for test_name in test_names:
		external_solver_avg_time = numpy.mean(results[test_name][solver_name])
		row['data'].append(round_to_2_digits(external_solver_avg_time))
	result_list.append(row)
# print result_list

# Replace strings in chart_template.html
f = open('chart_template.html')
new_f = open('chart.html', 'w')
for line in f:
	if '###TEST_NAMES###' in line:
		line = line.replace('###TEST_NAMES###', str(test_names))
	if '###DATA###' in line:
		line = line.replace('###DATA###', str(result_list))
	new_f.write(line)
f.close()
new_f.close()