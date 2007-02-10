"""
Implements a few data structures in SageX.

Written by Tom Boothby, 2007.  Free for any use.
"""
include '../ext/stdsage.pxi'
include '../ext/python.pxi'





cdef struct binary_tree_node
cdef struct binary_tree_node:
    int key
    binary_tree_node *left, *right
    void *value

cdef binary_tree_node *BinaryTreeNode(int key, object value):
    cdef binary_tree_node *t
    t = <binary_tree_node *>sage_malloc(sizeof(binary_tree_node))
    t.key = key
    t.left = NULL
    t.right = NULL
    Py_INCREF(value)
    t.value = <void *>value
    return t

cdef void free_binary_tree_node(binary_tree_node *self):
    if self.value != NULL:
        Py_DECREF(<object>self.value)
    sage_free(self)

cdef void binary_tree_dealloc(binary_tree_node *self):
    if self.left != NULL:
        binary_tree_dealloc(self.left)
        free_binary_tree_node(self.left)
    if self.right != NULL:
        binary_tree_dealloc(self.right)
        free_binary_tree_node(self.right)


cdef void binary_tree_insert(binary_tree_node *self, int key, object value):
    if self.key == key:
        return
    elif self.key > key:
        if self.left == NULL:
            self.left = BinaryTreeNode(key, value)
        else:
            binary_tree_insert(self.left, key, value)
    else:
        if self.right == NULL:
            self.right = BinaryTreeNode(key, value)
        else:
            binary_tree_insert(self.right, key, value)


cdef object binary_tree_get(binary_tree_node *self, int key):
    if self.key == key:
        return <object>self.value
    elif self.key > key:
        if self.left == NULL:
            return None
        else:
            return binary_tree_get(self.left, key)
    else:
        if self.right == NULL:
            return None
        else:
            return binary_tree_get(self.right, key)

cdef object binary_tree_delete(binary_tree_node *self, int key):
    cdef object t
    cdef binary_tree_node *cur
    if self.key > key:
        if self.left == NULL:
            return None
        elif self.left.key == key:
            t = <object>self.left.value
            self.left = binary_tree_left_excise(self.left)
            return t
        else:
            return binary_tree_delete(self.left, key)
    else:
        if self.right == NULL:
            return None
        elif self.right.key == key:
            t = <object>self.right.value
            self.right = binary_tree_right_excise(self.right)
            return t
        else:
            return binary_tree_delete(self.right, key)



cdef binary_tree_node *binary_tree_left_excise(binary_tree_node *self):
    cdef binary_tree_node *left, *cur
    if self.left == NULL:
        left = self.right
    elif self.right == NULL:
        left = self.left
    else:
        left = self.left
        cur = self.left
        while cur.right != NULL:
            cur = cur.right
        cur.right = self.left.right
    free_binary_tree_node(self)
    return left



cdef binary_tree_node *binary_tree_right_excise(binary_tree_node *self):
    cdef binary_tree_node *right, *cur
    if self.right == NULL:
        right = self.left
    elif self.left == NULL:
        right = self.right
    else:
        right = self.right
        cur = self.right
        while cur.left != NULL:
            cur = cur.left
        cur.left = self.right.left
    free_binary_tree_node(self)
    return right


cdef binary_tree_node *binary_tree_head_excise(binary_tree_node *self):
    cdef binary_tree_node *cur
    cdef int right
    # We have a pointer we're about to free.  Chances are, we'll never
    # see this pointer again.  Thus, it's least signifigant bit is
    # "random" enough to resist bias.
    right = (<int>self)&1
    if self.right == NULL:
        return self.left
    if self.left == NULL:
        return self.right
    if right:
        #move right branch to left, return left
        cur = self.left
        while cur.right != NULL:
            cur = cur.right
        cur.right = self.right
        cur = self.left
    else:
        #move left branch to right, return right
        cur = self.right
        while cur.left != NULL:
            cur = cur.left
        cur.left = self.left
        cur = self.right
    free_binary_tree_node(self)
    return cur

cdef class BinaryTree:
    """
    A simple binary tree with integer keys.
    """
    cdef binary_tree_node *head
    def __init__(BinaryTree self):
        self.head = NULL
    def __dealloc__(BinaryTree self):
        if self.head != NULL:
            binary_tree_dealloc(self.head)
            sage_free(self.head)

    def insert(BinaryTree self, int key, object value):
        """
        Inserts a key-value pair into the BinaryTree.  Duplicate keys are ignored.
        Example:
            sage: t = BinaryTree()
            sage: t.insert(1,1)
            sage: t.insert(0,0)
            sage: t.insert(2,2)
            sage: t.insert(0,1)
            sage: t.get(0)
            0
        """
        if self.head is NULL:
            self.head = BinaryTreeNode(key, value)
        else:
            binary_tree_insert(self.head, key, value)
    def delete(BinaryTree self, int key):
        """
        Removes a the node corresponding to key, and returns the value
        associated with it.
        Example:
            sage: t = BinaryTree()
            sage: t.insert(3,3)
            sage: t.insert(1,1)
            sage: t.insert(2,2)
            sage: t.insert(0,0)
            sage: t.insert(5,5)
            sage: t.insert(6,6)
            sage: t.insert(4,4)
            sage: t.delete(0)
            0
            sage: t.delete(3)
            3
            sage: t.delete(5)
            5
            sage: t.delete(2)
            2
            sage: t.delete(6)
            6
            sage: t.delete(1)
            1
            sage: t.delete(0)
            sage: t.get_max()
            4
            sage: t.get_min()
            4
        """
        cdef object r
        if self.head == NULL:
            return None
        elif self.head.key == key:
            r = <object>self.head.value
            self.head = binary_tree_head_excise(self.head)
            return r
        else:
            return binary_tree_delete(self.head, key)
    def get(BinaryTree self, int key):
        """
        Returns the value associated with the key given.
        Example:
            sage: t = BinaryTree()
            sage: t.insert(0,Matrix([[0,0],[1,1]]))
            sage: t.insert(0,1)
            sage: t.get(0)
            [0 0]
            [1 1]
        """
        if self.head == NULL:
            return <object>NULL
        else:
            return binary_tree_get(self.head, key)
    def contains(BinaryTree self, int key):
        """
        Returns True if a node with the given key exists
        in the tree, and False otherwise.
        Example:
            sage: t = BinaryTree()
            sage: t.contains(1)
            False
            sage: t.insert(1,1)
            sage: t.contains(1)
            True
        """
        if self.head == NULL:
            return False
        else:
            if binary_tree_get(self.head, key) is not None:
                return True
            else:
                return False
    def get_max(BinaryTree self):
        """
        Returns the value of the node with the maximal key value.
        """
        cdef binary_tree_node *cur
        if self.head == NULL:
            return None
        cur = self.head
        while cur.right != NULL:
            cur = cur.right
        return <object>cur.value
    def get_min(BinaryTree self):
        """
        Returns the value of the node with the minimal key value.
        """
        cdef binary_tree_node *cur
        if self.head == NULL:
            return None
        cur = self.head
        while cur.left != NULL:
            cur = cur.left
        return <object>cur.value
    def pop_max(BinaryTree self):
        """
        Returns the value of the node with the maximal key value,
        and removes that node from the tree.

        Example:
            sage: t = BinaryTree()
            sage: t.insert(4,'e')
            sage: t.insert(2,'c')
            sage: t.insert(0,'a')
            sage: t.insert(1,'b')
            sage: t.insert(3,'d')
            sage: t.insert(5,'f')
            sage: while not t.is_empty():
            ...    print t.pop_max()
            f
            e
            d
            c
            b
            a
        """
        cdef binary_tree_node *cur
        cdef object max
        if self.head == NULL:
            return None
        if self.head.right == NULL:
            max = <object>self.head.value
            cur = self.head.left
            free_binary_tree_node(self.head)
            self.head = cur
            return max
        cur = self.head
        while cur.right.right != NULL:
            cur = cur.right
        max = <object>cur.right.value
        cur.right = binary_tree_right_excise(cur.right)
        return max
    def pop_min(BinaryTree self):
        """
        Returns the value of the node with the minimal key value,
        and removes that node from the tree.
        Example:
            sage: t = BinaryTree()
            sage: t.insert(4,'e')
            sage: t.insert(2,'c')
            sage: t.insert(0,'a')
            sage: t.insert(1,'b')
            sage: t.insert(3,'d')
            sage: t.insert(5,'f')
            sage: while not t.is_empty():
            ...    print t.pop_min()
            a
            b
            c
            d
            e
            f
        """
        cdef binary_tree_node *cur
        cdef object min
        if self.head == NULL:
            return None
        if self.head.left == NULL:
            min = <object>self.head.value
            cur = self.head.right
            free_binary_tree_node(self.head)
            self.head = cur
            return min
        cur = self.head
        while cur.left.left != NULL:
            cur = cur.left
        min = <object>cur.left.value
        cur.left = binary_tree_left_excise(cur.left)
        return min
    def is_empty(BinaryTree self):
        """
        Returns True if the tree has no nodes.
        Example:
            sage: t = BinaryTree()
            sage: t.is_empty()
            True
            sage: t.insert(0,0)
            sage: t.is_empty()
            False
        """
        if self.head == NULL:
            return True
        else:
            return False
