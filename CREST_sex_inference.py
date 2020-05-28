#!/usr/bin/env python3

from __future__ import division
import numpy as np
import time
import sys



###########################
# Create an argument parser
###########################

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', dest = 'ifile', metavar = 'input_file_name', 
    required = True, help = 'Name of the input file')
parser.add_argument('-o', '--output', dest = 'ofile', metavar = 'output_file_name', 
    required = True, help = 'Name of the output file')
parser.add_argument('-m', '--map', dest = 'mfile', metavar = 'map_file_name', 
    required = True, help = 'Name of the genetic map file (should be in .simmap format)')
parser.add_argument('-b', '--bim', dest = 'bfile', metavar = 'bim_file',
        required = True, help = 'A PLINK .bim containing the dataset-specific map (should contain 22 autosomes)')

parser.add_argument('-w', '--window', metavar = 'window_size', 
    type = int, default = 500, help = 'Window size in kilobases. Default: 500 kb')
parser.add_argument('-k', '--keep', dest = 'keepfile', metavar = 'keep_file_name', 
    default = None, help = 'Name of the file containing samples to keep. If not provided, all samples are kept.')

args = parser.parse_args()

input_file = args.ifile
output_file = args.ofile
map_file = args.mfile
bim = args.bfile
keep_file = args.keepfile

# multiply by 1000 bp / kb * 0.5
window_size = args.window * 500



########################
# Functions to read data
########################

def read_input(file_name, keep = None):
  '''
  This function will reads in an input file for testing and converts it to a CREST format vector of vectors
  It can take in IBIS-style input, provided the rows correspond to segments and the fields are: ('id1', 'id2', 'chr', 'start (bp)', 'stop (bp)')
  If a keep_file is provided, it will remove rows not in the keep_file before returning this structured array.
  '''

  def reorder_pairs(pair):
    '''
    Enforce consistent ordering
    '''
    id1, id2 = pair
    sample1 = id1
    sample2 = id2
    if id2 > id1:
      sample1 = id2
      sample2 = id1
    return [sample1, sample2]

  # Read in the keep file
  if keep is not None:
    k = {}
    with open(keep) as kfile:
      for line in kfile:
        lk = reorder_pairs(line.strip().split()[:2])
        k[tuple(lk)] = True

  # Read in the segment file
  f = []

  with open(file_name) as file:
    for line in file:
      l = line.strip().split()[:5]

      if int(l[4]) <= int(l[3]):
        print('Segment end position must be greater than start position')
        sys.exit(1)
        
      pair = reorder_pairs(l[:2])
      l_modified = pair + l[2:]

      # need tuples to create a structured array
      if keep is None or tuple(pair) in k:
        # IBD2 segments may "interrupt" IBD1; if the previous segment would exactly coincide with the current one, merge the two
        if len(f) > 0 and tuple(f[-1][:2]) == tuple(l_modified[:2]) and f[-1][2] == l_modified[2] and abs(int(l_modified[3]) - int(f[-1][4])) < 2:
          intermediate_list = list(f[-1])
          intermediate_list[4] = l_modified[4]
          f[-1] = tuple(intermediate_list)
        else:
          f.append(tuple(l_modified))
      
  # create the structured array with appropriate fields
  f_in = np.array(f, dtype=[('id1', 'U64'), ('id2', 'U64'), ('chr', 'U64'), ('start', 'i8'), ('stop', 'i8')])

  return f_in



def read_simmap(file_name):
  '''
  Read in a .simmap file as a dictionary of structured arrays with fields: ('chrom', 'pos', 'male (cM)', 'female (cM)')
  '''
  f = []
  with open(file_name) as file:
    for line in file:
      # skip header lines
      if line[0] == "#":
        continue

      l = line.strip().split()
      
      # need tuples to create a structured array
      f.append(tuple(l))
      
  # create the structured array with appropriate fields
  f_in = np.array(f, dtype = [('chrom', 'U64'), ('pos', 'i8'), ('male', 'f8'), ('female', 'f8')])

  # split the simmap into chromosome maps
  D = dict()
  for chrom in np.unique(f_in['chrom']):
    D[chrom] = f_in[f_in['chrom'] == chrom]

  return D



def compute_bim_ends(file_name):
  '''
  Read in a .bim file and extract the ends of the bim map as a dictionary
  .bim files have fields: ('chr', 'SNP ID', 'pos (M)', 'bp (pos)', 'allele 1', 'allele 2')
  '''
  f = []
  old_line = None

  with open(file_name) as file:
    for line in file:
      l = line.strip().split()
      
      # if there is no previous line
      if not old_line:
        old_line = l
        f.append(tuple(old_line))

      # if the previous line is on the same chromosome
      if old_line[0] == l[0]:
        old_line = l

      # otherwise the line marks the end of one chromosome and the beginning of another
      else:
        f.append(tuple(old_line))
        f.append(tuple(l))
        old_line = l

    f.append(tuple(old_line))

  f_in = np.array(f, dtype = [('chrom', 'U64'), ('snp', 'U64'), ('morgan', 'f8'), ('pos', 'i8'), ('allele1', 'U64'), ('allele2', 'U64')])
  
  D = dict()
  for chrom in np.unique(f_in['chrom']):
    D[chrom] = [f_in[f_in['chrom'] == chrom]['pos'][0], f_in[f_in['chrom'] == chrom]['pos'][-1]]

  return D



#######################################
# Functions for manipulating logarithms
#######################################

def log1p_10(x):
  '''
  Calculate log_10(1 + x) accurately, exists as a function only to beautify code
  '''
  return np.log1p(x) / np.log(10)



def log_sum(list_of_logs):
  '''
  Given a list of base 10 log values: [ log_10(x0), ..., log_10(xN) ], find the log-sum: log_10(x0 + ... + xN)
  '''
  a = list_of_logs[0]

  x = 1
  while x < len(list_of_logs):
    b = list_of_logs[x]
    if a > b:
      a = a + log1p_10(10 ** (b - a))
    else:
      a = b + log1p_10(10 ** (a - b))
    x += 1

  return a



def log_product(list_of_logs):
  '''
  Given a list of base 10 log values: [ log_10(x0), ..., log_10(xN) ], find the log-product: log_10(x0 * ... * xN)
  '''
  return sum(list_of_logs) 



def log_factorial(n):
  '''
  For an integer n, compute log_10(n!)
  '''
  fact_n = 0
  if n == 0 or n == 1:
    return fact_n

  for i in range(1, n + 1):
    fact_n += np.log10(i)
  return fact_n



def log_poisson(events, rate):
  '''
  Compute the log_10 of the Poisson p.m.f. for the probability of, e.g. k events with rate r
  f(k; r) = r^k * e^-r / k!
  log_10(f(k; r)) = k*log_10(r) - r*log_10(e) - log_factorial(k)
  '''
  if rate == 0:
    return 0
  elif rate < 0:
    print('Warning, Poisson rate negative. Check segment.')
  else:
    return (events * np.log10(rate)) - (rate * np.log10(np.e)) - log_factorial(events)



#####################################
# Functions to handle the genetic map
#####################################
def linear_interpolation(x, p0, p1):
  '''
  Perform linear interpolation to find cM values for points not attested on the genetic map
  Input is a point x not necessarily on the map, and p0, p1, two flanking points
  Output px, a point in the format of p0, p1, i.e.: ('chr', 'pos (bp)', 'cM (male)', 'cM (female)')
  '''
  chrom0, x0, male0, female0 = p0
  chrom1, x1, male1, female1 = p1
  
  if chrom0 != chrom1:
    print('Warning: chromosome mismatch, returning p0')
    return p0

  # don't interpolate if p0 and p1 overlap
  elif x0 == x1:
    return p0
  
  # don't interpolate if x falls on or outside the bounds for interpolation
  elif x <= x0:
    return p0
  
  elif x >= x1:
    return p1
 
  else:
    # get interpolation for both maps
    p0_cM = [male0, female0]
    p1_cM = [male1, female1]
  
    # will expand via a loop
    px = [chrom0, x]

    for i in range(2):
      # first male, then female
      y0 = p0_cM[i]
      y1 = p1_cM[i]

      # calculate y
      y = y0 + (x - x0)*(y1 - y0)/(x1 - x0)
      px.append(y)

    return px



def pos_search(pos, simmap):
  '''
  Perform a quasi-binary search for a position 'pos' in a 'sex' sex-specific single-chromosome genetic map 'simmap'
  In this context simmap is an individual structured array, which is an entry in the dictionary produced by read_simmap()
  Output: interpolated position of the same format ('chr', 'pos (bp)', 'cM (male)', 'cM (female)') OR the index of the closest position on the input map
  Credit Amy Williams, adapted to python from C by me
  '''
  # get the appropriate chromosome and sex of the genetic map
  gmap = simmap['pos']
  
  # starting position for the search
  map_index = 0
  
  # index bounds within which the position falls
  left, right = [0, len(gmap)]
  
  if pos >= gmap[right - 1]:
    map_index = right - 1
  
  elif pos <= gmap[left]: 
    map_index = left
  
  else:
    while True:
      if right - left == 1:
        map_index = left
        break
      
      mid = (left + right) // 2
      if gmap[mid] < pos:
        left = mid
      elif gmap[mid] > pos:
        right = mid
      else:
        # equal: exact map position
        map_index = mid;
        right = mid
        left = mid
        break;
  
  if right == len(gmap):
    right -= 1
  
  # perform linear interpolation between left and right
  p0 = simmap[left]
  p1 = simmap[right]
  return linear_interpolation(pos, p0, p1)



def get_windows(bim, simmap, segment):
  '''
  Find windows overlapping an IBD segment on the genetic map and get their lengths in Morgans
  The end_processing input dictates how segments approaching the map ends should be handled
  '''
  # get the segment information
  chrom, start, stop = segment
 
  # list of positions to pass (suffix 0 means the head, 1 means the tail)
  start_0 = start - window_size
  start_1 = start + window_size
  stop_0 = stop - window_size
  stop_1 = stop + window_size
  
  # genetic map
  gmap = simmap[chrom]
  gmap_ends = (gmap['pos'][0], gmap['pos'][-1])
  bim_ends = (bim[chrom][0], bim[chrom][1])
  usable_gmap = gmap
  usable_gmap_ends = (max(bim_ends[0], gmap_ends[0]), min(bim_ends[1], gmap_ends[1]))

  # forbid windows from exceeding or reaching the ends of the map
  if start <= usable_gmap_ends[0]:
    start_0 = usable_gmap_ends[0]
    start_1 = usable_gmap_ends[0]
  
  if stop >= usable_gmap_ends[1]:
    stop_0 = usable_gmap_ends[1]
    stop_1 = usable_gmap_ends[1]

  # the interval is created from the two windows  
  interval_0 = start_1
  interval_1 = stop_0
  
  pos_list = [(start_0, start_1), (interval_0, interval_1), (stop_0, stop_1)]

  lengths = dict()
  lengths['male'] = []
  lengths['female'] = []

  for (x, y) in pos_list:
    if y < x:
      return {'male':[0, 0, 0], 'female':[0, 0, 0]}
    x_chrom, x_pos, x_cM_male, x_cM_female = pos_search(x, usable_gmap)
    y_chrom, y_pos, y_cM_male, y_cM_female = pos_search(y, usable_gmap)
    if y_pos < x_pos:
      return {'male':[0, 0, 0], 'female':[0, 0, 0]}
    lengths['male'].append((y_cM_male - x_cM_male)/100)
    lengths['female'].append((y_cM_female - x_cM_female)/100)

  return lengths



def segment_gaps(bim, simmap, all_segments):
  '''
  For all segments shared, get the lengths of the non-segment gaps
  '''

  # for every chromosome, get the length of every gap on both sex specific maps
  segs = dict()
  gaps = dict()
  gaps['male'] = []
  gaps['female'] = []

  for chrom in simmap:
    # get the segments sorted by start position and separated by chromosome
    segs[chrom] = [(seg[1], seg[2]) for seg in all_segments if seg[0] == chrom]
    segs[chrom].sort()

    # define the most conservative map edges
    gap_list = []
    bim_ends = (bim[chrom][0], bim[chrom][1])
    map_ends = (simmap[chrom]['pos'][0], simmap[chrom]['pos'][-1])
    gmap_ends = (max(bim_ends[0], map_ends[0]), min(bim_ends[1], map_ends[1]))

    # gaps must have positive lengths
    # the first gap falls between the first segment and the map start
    gap_start = gmap_ends[0]

    # if there are no segments on the chromosome, the gap is exactly the map ends
    if segs[chrom] == []:
      gap_end = gmap_ends[1]
      gap_list.append((gap_start, gap_end))

    else:
      gap_end = segs[chrom][0][0] - window_size
      if gap_end > gap_start:
        gap_list.append((gap_start, gap_end))

      # fill in between each segment and the next
      for s in range(len(segs[chrom]) - 1):
        gap_start = segs[chrom][s][1] + window_size
        gap_end = segs[chrom][s + 1][0] - window_size
        if gap_end > gap_start:
          gap_list.append((gap_start, gap_end))

      # the last gap falls between the last segment and the map end
      gap_start = segs[chrom][-1][1] + window_size
      gap_end = gmap_ends[1]
      if gap_end > gap_start:
        gap_list.append((gap_start, gap_end))

    # x and y are the boundaries of the gap (not to be confused with the sex chromosomes)
    for (x, y) in gap_list:
      x_chrom, x_pos, x_cM_male, x_cM_female = pos_search(x, simmap[chrom])
      y_chrom, y_pos, y_cM_male, y_cM_female = pos_search(y, simmap[chrom])
      if x_pos < y_pos:
        gaps['male'].append((y_cM_male - x_cM_male)/100)
        gaps['female'].append((y_cM_female - x_cM_female)/100)

  return gaps



#################################
# Functions to compute LOD scores
#################################

def logp_k_crossovers(k, segment):
  '''
  Compute the log_10 probability of observing k crossovers in a segment
  The segment should be in the format of an entry in the dictionary returned by get_windows(...)
  '''
  logp_start = log_poisson(1, segment[0])
  logp_interval = log_poisson(k, segment[1])
  logp_stop = log_poisson(1, segment[2])

  return log_product([logp_start, logp_interval, logp_stop])



def logp_gaps(gaps):
  '''
  Get the probability of no internal recombinations for each gap
  '''
  l_P_gaps = [log_poisson(0, gap) for gap in gaps]
  return l_P_gaps



# These two functions calculate LODs
# input is the dictionary produced in get_LODs: pairs[('id1', 'id2')] = [ seg1:(chr, start, stop), ..., segN:(chr, start, stop)]
# output is a pair of LOD scores (female_LOD, male_LOD)

def LOD_gp(bim, simmap, data, key):
  
  list_logp_female = []
  list_logp_male = []

  for segment in data:
    # get the window lengths
    lengths = get_windows(bim, simmap, segment)
    list_logp_female.append(logp_k_crossovers(0, lengths['female']))
    list_logp_male.append(logp_k_crossovers(0, lengths['male']))

  all_gaps = segment_gaps(bim, simmap, data)
  list_logp_gaps_female = logp_gaps(all_gaps['female'])
  list_logp_gaps_male = logp_gaps(all_gaps['male'])

  logp_female = log_product(list_logp_female + list_logp_gaps_female)
  logp_male = log_product(list_logp_male + list_logp_gaps_male)

  LOD = logp_female - logp_male
  return LOD



def LOD_hs(bim, simmap, data, key):

  list_logp_female = []
  list_logp_male = []

  for segment in data:
    # get the window lengths
    lengths = get_windows(bim, simmap, segment)
    # adjust effective lengths for the 2 HS meioses
    lengths['female'] = [2 * x for x in lengths['female']]
    lengths['male'] = [2 * x for x in lengths['male']]

    list_logp_female.append(logp_k_crossovers(0, lengths['female']))
    list_logp_male.append(logp_k_crossovers(0, lengths['male']))

  all_gaps = segment_gaps(bim, simmap, data)
  all_gaps['female'] = [2 * x for x in all_gaps['female']]
  all_gaps['male'] = [2 * x for x in all_gaps['male']]

  list_logp_gaps_female = logp_gaps(all_gaps['female'])
  list_logp_gaps_male = logp_gaps(all_gaps['male'])

  logp_female = log_product(list_logp_female + list_logp_gaps_female)
  logp_male = log_product(list_logp_male + list_logp_gaps_male)

  LOD = logp_female - logp_male
  return LOD



def get_LODs(pairs, bim, simmap):
  '''
  For a dictionary with format:
  
  pairs[('id1', 'id2')] = [ seg1:(chr, start, stop),
                            seg2:(chr, start, stop),
                            ...
                            segN:(chr, start, stop) ]

  compute LOD scores for each pair for each of the three relationship types
  '''
  # the dictionary LODs holds the segment numbers (line count) and the LOD scores
  LODs = {}

  for pair in pairs:
    LODs[pair] = np.array((0, 0, 0), dtype = [('segnum', 'i8'), ('gpLOD', 'f8'), ('hsLOD', 'f8')])
    data = pairs[pair]
    segnum = len(data)
    LODs[pair]['segnum'] = segnum
    LODs[pair]['gpLOD'] = LOD_gp(bim, simmap, data, pair)
    LODs[pair]['hsLOD'] = LOD_hs(bim, simmap, data, pair)
  return LODs



##########################
# Function to write output
##########################

def write_output(input_struct, file_name, bim, simmap):
  '''
  Write out results of IBD inference
  '''

  # construct a dictionary of unique pairs
  pairs_uniq = np.unique(input_struct[['id1', 'id2']])
  pairs_dict = dict()
  
  for x in range(len(pairs_uniq)):
    pair = tuple(pairs_uniq[x])
    pairs_dict[pair] = []
  
  for y in range(len(input_struct)):
    pair = tuple(input_struct[['id1', 'id2']][y])
    data = input_struct[['chr', 'start', 'stop']][y]
    pairs_dict[pair].append(data)

  # compute LOD scores for all pairs
  out = get_LODs(pairs_dict, bim, simmap)
  
  # write output
  with open(file_name, 'w') as file:
    for pair in out:
      line_to_write = list(pair)
      data_to_add = [str(i) for i in out[pair].tolist()]
      line_to_write += data_to_add
      file.write('\t'.join(line_to_write) + '\n')
  return out



###################################
# Read input, get map, write output
###################################

def main():
  t0 = time.time()
  in_file = read_input(input_file, keep_file)
  bim_ends = compute_bim_ends(bim)
  simmap = read_simmap(map_file)
  write_output(in_file, output_file, bim_ends, simmap)
  t1 = time.time()
  print('Total run-time is:', t1 - t0)


if __name__ == '__main__':
  main()
