from random import choice
import time
import memory_profiler
import matplotlib.pyplot as plt


def weightedchoice(items): # this doesn't require the numbers to add up to 100
    return choice("".join(x * y for x, y in items))

def String(length):

   DNA=""
   for count in range(length):
      DNA+=weightedchoice([("C", 10), ("G", 20), ("A", 40), ("T", 30)])
   return DNA

# Python program for KMP Algorithm
def KMPSearch(pat, txt):
	M = len(pat)
	N = len(txt)

	# create lps[] that will hold the longest prefix suffix
	# values for pattern
	lps = [0]*M
	j = 0 # index for pat[]

	# Preprocess the pattern (calculate lps[] array)
	computeLPSArray(pat, M, lps)

	i = 0 # index for txt[]
	while i < N:
		if pat[j] == txt[i]:
			i += 1
			j += 1

		if j == M:
			print ("Found pattern at index {} by KMP Search".format(str(i-j)))
			j = lps[j-1]

		# mismatch after j matches
		elif i < N and pat[j] != txt[i]:
			# Do not match lps[0..lps[j-1]] characters,
			# they will match anyway
			if j != 0:
				j = lps[j-1]
			else:
				i += 1

def computeLPSArray(pat, M, lps):
	len = 0 # length of the previous longest prefix suffix

	lps[0] # lps[0] is always 0
	i = 1

	# the loop calculates lps[i] for i = 1 to M-1
	while i < M:
		if pat[i]== pat[len]:
			len += 1
			lps[i] = len
			i += 1
		else:
			# This is tricky. Consider the example.
			# AAACAAAA and i = 7. The idea is similar
			# to search step.
			if len != 0:
				len = lps[len-1]

				# Also, note that we do not increment i here
			else:
				lps[i] = 0
				i += 1

NO_OF_CHARS = 256
 
def badCharHeuristic(string, size):
    '''
    The preprocessing function for
    Boyer Moore's bad character heuristic
    '''
 
    # Initialize all occurrence as -1
    badChar = [-1]*NO_OF_CHARS
 
    # Fill the actual value of last occurrence
    for i in range(size):
        badChar[ord(string[i])] = i;
 
    # return initialized list
    return badChar
 
def Boyer_Moore_search(txt, pat):
    '''
    A pattern searching function that uses Bad Character
    Heuristic of Boyer Moore Algorithm
    '''
    m = len(pat)
    n = len(txt)
 
    # create the bad character list by calling
    # the preprocessing function badCharHeuristic()
    # for given pattern
    badChar = badCharHeuristic(pat, m)
 
    # s is shift of the pattern with respect to text
    s = 0
    while(s <= n-m):
        j = m-1
 
        # Keep reducing index j of pattern while
        # characters of pattern and text are matching
        # at this shift s
        while j>=0 and pat[j] == txt[s+j]:
            j -= 1
 
        # If the pattern is present at current shift,
        # then index j will become -1 after the above loop
        if j<0:
            print("Found pattern at index {} by Boyer_Moore_searchˇ".format(s))
 
           
            s += (m-badChar[ord(txt[s+m])] if s+m<n else 1)
        else:
           
            s += max(1, j-badChar[ord(txt[s+j])])


def Naive_search(pat, txt):
    M = len(pat)
    N = len(txt)
 
    # A loop to slide pat[] one by one */
    for i in range(N - M + 1):
        j = 0
 
        # For current index i, check
        # for pattern match */
        while(j < M):
            if (txt[i + j] != pat[j]):
                break
            j += 1

        if(j == M):
            print("Pattern found at index {} by Naive search algorithm".format(i))




def main():
    KMPSearch_memory = []
    Boyer_Moore_search_memory = []
    Naive_search_memory = []
    indexs = []
    
    KMPSearch_time = []
    Boyer_Moore_search_time = []
    Naive_search_time = []
    
    print('test for correctness of algorithm : ')
    print()
    txt = String(20)
    pat = String(3)
    print("this is random generated gene : {}".format(txt))
    print("this is part of gene we are looking for : {}".format(pat))
    KMPSearch(pat, txt)
    Boyer_Moore_search(txt, pat)
    Naive_search(pat, txt)
    
    print()
    print('------------------------------------------')

    KMPSearch(pat, txt)
    Boyer_Moore_search(txt, pat)
    Naive_search(pat, txt)
    indexs = []
    txt_length = [10,100,1000,10000]
    pat_length = [2,4,6,8]
    
    for i in range(len(txt_length)):
        
        for j in range(len(txt_length)):
            index = []
            index.append(txt_length[i])
            index.append(pat_length[j])
            indexs.append(index)
            txt = String(txt_length[i])
            pat = String(pat_length[i])



            
            init = memory_profiler.memory_usage()
            start = time.time()
            
            KMPSearch(pat, txt)
            
            end = time.time()
            finish = memory_profiler.memory_usage()
            execution_time = end - start 
            space_usage = finish[0] - init[0]
            KMPSearch_memory.append(space_usage)
            KMPSearch_time.append(execution_time)
        
        # -------------------------- -----------------------------    
            init = memory_profiler.memory_usage()
            start = time.time()
            
            Boyer_Moore_search(txt, pat)
            
            end = time.time()
            finish = memory_profiler.memory_usage()
            execution_time = end - start 
            space_usage = finish[0] - init [0]
            Boyer_Moore_search_memory.append(space_usage)
            Boyer_Moore_search_time.append(execution_time)
        
        # -------------------------- -----------------------------   
            init = memory_profiler.memory_usage()
            start = time.time()
            
            Naive_search(pat, txt)
            
            end = time.time()
            finish = memory_profiler.memory_usage()
            execution_time = end - start 
            space_usage = finish[0] - init [0]
            Naive_search_memory.append(space_usage)
            Naive_search_time.append(execution_time)
    number = len(KMPSearch_memory)    
    index = range(0 , number)
    plt.figure(figsize=(15, 6))
    plt.subplot(2,1,1)

    plt.xticks(range(0,len(indexs)), indexs)
    plt.plot(index , KMPSearch_memory , label = 'KMPSearch_memory')
    plt.plot(index , Boyer_Moore_search_memory , label = 'Boyer_Moore_search_memory')
    plt.plot(index , Naive_search_memory , label = 'Naive_search_memory')
    plt.legend()
    
    plt.subplot(2,1,2)

    plt.xticks(range(0,len(indexs)), indexs)
    plt.plot(index , KMPSearch_time , label = 'KMPSearch_time')
    plt.plot(index , Boyer_Moore_search_time , label = 'Boyer_Moore_search_time')
    plt.plot(index , Naive_search_time , label = 'Naive_search_time')
    
    plt.legend()
  
    
    
    print()
    print(KMPSearch_time)
    
if __name__ == '__main__':
    main()
    
    
    