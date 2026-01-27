# Lab Assignment 3: Practicing Operations.

## Imports

import time
import random

## 1. Practicing operations on a Python list: append(), clear(), copy(), count(), extend(), index(), insert(), pop(), remove(), reverse(), sort().

def list_operations():
    my_list = [random.randint(0, 1000) for _ in range(10)] # Initial random list of 10 elements
    print("Original list:", my_list)

    # Append
    my_list.append(120)
    print("After append(120):", my_list)

    # Count the occurrence of a specific element, guessing 5 here
    count_5 = my_list.count(5)
    print("Count of element 5:", count_5)

    # Extend
    my_list.extend([12, 13]) # Extending with another list
    print("After extend([12, 13]):", my_list)

    # Index
    index_of_3 = my_list.index(3) if 3 in my_list else print("3 is not in list to find index.")
    print("Index of 3:", index_of_3)

    # Insert
    my_list.insert(0, 0) # Insert 0 at index 0
    print("After insert(0, 0):", my_list)

    # Pop
    popped_value = my_list.pop() 
    print("After pop():", my_list, "\nPopped value:", popped_value)

    # Remove -  Guarded to avoid error if 5 not in list
    if 5 in my_list:
        my_list.remove(5)
        print("After remove(5):", my_list)
    else:
        print("5 not in list to remove.")

    # Reverse
    my_list.reverse() # Reverse the list
    print("After reverse():", my_list)

    # Sort
    my_list.sort() # Sort the list
    print("After sort():", my_list)

    # Copy
    copied_list = my_list.copy()
    print("Copied list:", copied_list)

    # Clear
    my_list.clear() # We should see an empty list after this
    print("After clear():", my_list) 

## 2. Generate a random list L of length N. Compare the performance of the following operations by measuring the times with increasing N
## N: N, 2N, 4N, 8N,... explain what you observe.

def perf_comparison():
    """
    Starting with initial N of 100, then doubling it with each iteration up to a maximum of 10N.
    The operations performed are:
    1. Cost of removing elements from lists: L.pop() vs L.pop(0)
    2. Slicing vs explicit copy of a list: A = L[:] vs B = list(L)
    3. In-place vs out-of-place modification: L.reverse() vs R = L[::-1]
    """
    initial_N = 100
    # max_multiplier = 10
    operation_times = {
        "pop_end": [],
        "pop_start": [],
        "slice_copy": [],
        "list_copy": [],
        "inplace_reverse": [],
        "outofplace_reverse": []
    } # A dictionary to hold times for each operation

    #for multiplier in range(1, max_multiplier + 1):
    for multiplier in [1, 2, 4, 8, 10]:    
        N = initial_N * multiplier
        L = [random.randint(1, 1000) for _ in range(N)]

        # Measure pop() from end
        start_time = time.time()
        for _ in range(N):
            L.pop()
        operation_times["pop_end"].append(time.time() - start_time)

        # Recreate list
        L = [random.randint(1, 1000) for _ in range(N)]

        # Measure pop(0) from start
        start_time = time.time()
        for _ in range(N):
            L.pop(0)
        operation_times["pop_start"].append(time.time() - start_time)

        # Recreate list
        L = [random.randint(1, 1000) for _ in range(N)]

        # Measure slicing copy
        start_time = time.time()
        A = L[:]
        operation_times["slice_copy"].append(time.time() - start_time)

        # Measure list() copy
        start_time = time.time()
        B = list(L)
        operation_times["list_copy"].append(time.time() - start_time)

        # Measure in-place reverse
        start_time = time.time()
        L.reverse()
        operation_times["inplace_reverse"].append(time.time() - start_time)

        # Recreate list
        L = [random.randint(1, 1000) for _ in range(N)]

        # Measure out-of-place reverse
        start_time = time.time()
        R = L[::-1]
        operation_times["outofplace_reverse"].append(time.time() - start_time)

    
    # Print results
    for operation, times in operation_times.items():
        print(f"{operation}: {times}")

    # Checking the changes in time with increaasing N in each list
    # This involves checking the ratio of times between consecutive N values for the operations in operation_times
    print("\nTime Ratios for Increasing N:")
    # for all elements in operation_times
    for operation, times in operation_times.items():
        print(f"\nOperation: {operation}")
        for i in range(1, len(times)):
            ratio = times[i] / times[i - 1] if times[i - 1] != 0 else float('inf')
            print(f"  Ratio of time for N={initial_N * (2**i)} to N={initial_N * (2**(i-1))}: {ratio:.2f}")
       
 

## Execution
if __name__ == "__main__":
    print("List Operations:")
    list_operations()
    print("\nPerformance Comparison:")
    perf_comparison()
