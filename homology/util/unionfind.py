'''
unionfind.py

http://code.activestate.com/recipes/215912/

A class that implements the Union Find data structure and algorithm.  This
data structure allows one to find out which set an object belongs to, as well
as join two sets.

The algorithm's performance, given m union/find operations of any ordering, on
n elements has been shown to take log* time per operation, where log* is
pronounced log-star, and is the INVERSE of what is known as the Ackerman
function, which is given below:
A(0) = 1
A(n) = 2**A(n-1)

I include the functions to be complete.  Note that we can be 'inefficient'
when performing the inverse ackerman function, as it will only take a maximum
of 6 iterations to perform; A(5) is 65536 binary digits long (a 1 with 65535
zeroes following).  A(6) is 2**65536 binary digits long, and cannot be
represented by the memory of the entire universe.


The Union Find data structure is not a universal set implementation, but can
tell you if two objects are in the same set, in different sets, or you can
combine two sets.
ufset.find(obja) == ufset.find(objb)
ufset.find(obja) != ufset.find(objb)
ufset.union(obja, objb)


This algorithm and data structure are primarily used for Kruskal's Minimum
Spanning Tree algorithm for graphs, but other uses have been found.

August 12, 2003 Josiah Carlson
Modified: Aug. 2017, Jose Licon
'''

import numpy as np


class UnionFind:

    def __init__(self):
        '''
        Create an empty union find data structure.
        '''
        self.num_weights = {}
        self.parent_pointers = {}
        self.num_to_objects = {}
        self.objects_to_num = {}
        self.__repr__ = self.__str__

    def insert_objects(self, objects, **kwargs):
        '''
        Insert a sequence of objects into the structure.  All must be Python hashable.
        '''
        for object in objects:
            self.find(object)

    def find(self, object, **kwargs):
        '''
        Find the root of the set that an object is in.
        If the object was not known, will make it known, and it becomes its own set.
        Object must be Python hashable.
        '''
        value = 1
        
        if object not in self.objects_to_num:
            print "new element:" , object
            obj_num = len(self.objects_to_num)
            self.num_weights[obj_num] = value
            self.objects_to_num[object] = obj_num
            self.num_to_objects[obj_num] = object
            self.parent_pointers[obj_num] = obj_num
            print "returning: ", object, " of ", type(object)
            return object
        stk = [self.objects_to_num[object]]
        par = self.parent_pointers[stk[-1]] 
        while par != stk[-1]:
            stk.append(par)
            par = self.parent_pointers[par]
        for i in stk:
            self.parent_pointers[i] = par
        return self.num_to_objects[par]

    def union(self, object1, object2, elder_rule=False):
        '''
        Combine the sets that contain the two objects given.
        Both objects must be Python hashable.
        If either or both objects are unknown, will make them known, and combine them.
        '''
        o1p = self.find(object1)
        o2p = self.find(object2)
        if o1p != o2p:
            on1 = self.objects_to_num[o1p]
            on2 = self.objects_to_num[o2p]
            if elder_rule:
                w1 = self.num_weights[on1]
                w2 = self.num_weights[on2]
                if on2 < on1:
                    o1p, o2p, on1, on2, w1, w2 = o2p, o1p, on2, on1, w2, w1
            else:
                w1 = self.num_weights[on1]
                w2 = self.num_weights[on2]
                if w1 < w2:
                    o1p, o2p, on1, on2, w1, w2 = o2p, o1p, on2, on1, w2, w1
            self.num_weights[on1] = w1 + w2
            del self.num_weights[on2]
            self.parent_pointers[on2] = on1

    def __str__(self):
        '''
        Included for testing purposes only.
        All information needed from the union find data structure can be attained
        using find.
        '''
        sets = {}
        for i in xrange(len(self.objects_to_num)):
            sets[i] = []
        for i in self.objects_to_num:
            sets[self.objects_to_num[self.find(i)]].append(i)
        out = []
        for i in sets.itervalues():
            if i:
                out.append(repr(i))
        return ', '.join(out)


class UnionFindOpt(object):

    def __init__(self, n):
        """
        n : int
        """
        self.parent = np.array(range(n), dtype=np.uint32)
        self.size   = np.ones(n, dtype=np.uint32)
        # Initialize using some value
        self.keys = {k : 0 for k in range(n)}

    def find(self, i):
        """
        i : int
        return p : int
        """
        a = i
        b = self.parent[i]
        if a == b:
            # Node is its own parent, do nothing
            return i
        # Go up the hierarchy until root is found
        while a != b:
            a = b
            b = self.parent[a]
        # Root is now stored in b
        p = b
        # Return to starting level, go back up pointing everything to the root along the way.
        a = i
        b = self.parent[i]
        while a != b:
            self.parent[a]  = p
            a = b
            b = self.parent[a]
        # Return root
        return p

    def union(self, i, j):
        """
        i : int
        j : int
        return : void
        """
        a = self.find(i)
        b = self.find(j)
        if a == b:
            return
        if self.size[a] > self.size[b]:
            self.parent[b] = a
            self.size[a] += self.size[b]
            self.size[b] = 0
        else:
            self.parent[a] = b
            self.size[b] += self.size[a]
            self.size[a] = 0


if __name__ == '__main__':

    uf = UnionFindOpt(11)
    print("parents at start: %s" % uf.parent)
    uf.union(0, 1)
    uf.union(2, 1)
    uf.union(3, 1)
    uf.union(4, 3)
    uf.union(6, 5)
    uf.union(7, 5)
    uf.union(10, 5)
    uf.union(9, 7)
    uf.union(8, 7)

    # p = np.array([ 3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3, 13,  3,  3,  3])
    # assert(all(uf.parent == p))
    r = np.array([1, 5, 1, 1, 1, 6, 1, 1, 1, 1, 1])
    p = np.array([1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5])
    assert(all(p == uf.parent))
    assert(all(r == uf.size))
