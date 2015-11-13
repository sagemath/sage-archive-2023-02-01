
## Dinakar Muthiah June 3, 2014

""" Sample usage
W = CoxeterGroup(['A',5])
load('BraidMoveCalculator.sage'); B = BraidMoveCalculator(W)
view(B.chain_of_reduced_words((1,2,1,3,2,1,4,3,2,1,5,4,3,2,1),(5,4,5,3,4,5,2,3,4,5,1,2,3,4,5)))

"""
from sage.combinat.root_system.cartan_type import CartanType


class BraidMoveCalculator(object):
    def __init__(self,coxeter_group):
        self.coxeter_group = coxeter_group
        self.simple_reflections = self.coxeter_group.simple_reflections()
           
    def multiply_word(self,word):
        return self.coxeter_group.from_reduced_word(word)
        #return prod([self.s[i] for i in word])
    
    def is_reduced(self,word):
        return (len(word) == self.multiply_word(word).length())

    def braid_word(self,i,j):
        s = self.simple_reflections
        # m_ij = (s[i]*s[j]).order()  #For some reason this is no longer implemented
        m_ij = self.coxeter_group.coxeter_matrix()[i-1,j-1]
        num_pairs = int(m_ij)/int(2)
        extra_term = int(m_ij)%2
        return (i,j)*num_pairs + (i,)*extra_term

    def _apply_put_in_front_recur_step(self, k, input_word, coxeter_matrix_entry):        
        #print "k:", k
        #print "input_word:", input_word
        #print "coxeter_matrix_entry:", coxeter_matrix_entry
        i = input_word[0]
        def partial_braid_word(length,swap=False,i=i,k=k):
            if swap:
                i,k = k,i
            running_braid_word = []
            for counter in xrange(length):
                if counter % 2 == 0:
                    running_braid_word.append(i)
                else:
                    running_braid_word.append(k)
            return tuple(running_braid_word)
            
        current_last_word = input_word
        current_first_letter = k
        output_word_list = [current_last_word] 
        for counter in xrange(1, coxeter_matrix_entry):
            current_word_list = self.put_in_front(current_first_letter, current_last_word[1:])
            output_word_list += [partial_braid_word(counter) + word 
                                 for word in current_word_list[1:]]
            if current_first_letter == k:
                current_first_letter = i
            else:
                current_first_letter = k
            current_last_word = current_word_list[-1]
            #print "counter:", counter
            #print "current_last_word", current_last_word
            #print "current_first_letter", current_first_letter
        if i != k:
            output_word_list += [partial_braid_word(coxeter_matrix_entry,swap=True) +
                                 current_last_word[1:]]
        return output_word_list

    def put_in_front(self, k, input_word):
        """
        Return a list of reduced words beginning with 
        ``input_word`` and ending with a reduced word whose first letter 
        is ``k``. There still remains an issue with 0 indices.
        """
        i = input_word[0]
        if i == 0 or k == 0:
            raise NotImplementedError
        entry = self.coxeter_group.coxeter_matrix()[i, k]
        return self._apply_put_in_front_recur_step(k, input_word, entry)

    def chain_of_reduced_words(self, start_word, end_word):
        if start_word == end_word:
            return [start_word]
        k = end_word[0]
        first_word_list = self.put_in_front(k, start_word)
        first_last_word = first_word_list[-1]
        return (first_word_list[:-1] + 
                [ (k,) + word for word in
                  self.chain_of_reduced_words(first_last_word[1:],
                                                   end_word[1:])])


# Maybe we need to be more specific, and pass not the Cartan type, but the root lattice?
# TODO: Move to PBW_data?
def enhance_braid_move_chain(braid_move_chain, cartan_type):
    r"""
    Return a list of tuples that records the data of the long words in
    ``braid_move_chain`` plus the data of the intervals where the braid moves
    occur and the data of the off-diagonal entries of the `2 \times 2` Cartan
    submatrices of each braid move.

    INPUT:

    - ``braid_move_chain`` -- a chain of reduced words in the Weyl group
      of ``cartan_type``
    - ``cartan_type`` -- a finite-type CartanType

    OUTPUT:
    A list of 3-tuples ``(reduced_word,interval_of_change,cartan_sub_matrix)``
    where 
    - ``reduced_word`` agrees with the elements of ``braid_move_chain``
    - ``interval_of_change`` is the (half-open) interval of indices where the 
      braid move occurs. This is `None` for the first tuple.
    - ``cartan_sub_matrix`` is the off-diagonal entries of the `2 \times 2`
      submatrix of the cartan matrix corresponding to the braid move.
      This is `None` for the first tuple.

    For a matrix::

        [2 a]
        [b 2]

    the ``cartan_sub_matrix`` is the pair ``(a, b)``.
    """
    cartan_type = CartanType(cartan_type)
    if not cartan_type.is_finite():
        raise ValueError("cartan_type must be finite type")
    cartan_matrix = cartan_type.cartan_matrix()
    output_list = []
    output_list.append( (braid_move_chain[0],None,None) )
    for k,current_word in enumerate(braid_move_chain[1:]):
        previous_word = braid_move_chain[k]
        interval_of_change = diff_interval(previous_word,current_word)
        i = previous_word[interval_of_change[0]] - 1  # -1 for indexing
        j = current_word[interval_of_change[0]] - 1  # -1 for indexing
        cartan_sub_matrix = (cartan_matrix[i,j], cartan_matrix[j,i])
        output_list.append( (current_word,interval_of_change,cartan_sub_matrix) )
    return output_list


def enhanced_braid_chain_test():
    braid_chain = [(1, 2, 1, 3, 2, 1),
                   (1, 2, 3, 1, 2, 1),
                   (1, 2, 3, 2, 1, 2),
                   (1, 3, 2, 3, 1, 2),
                   (3, 1, 2, 3, 1, 2),
                   (3, 1, 2, 1, 3, 2),
                   (3, 2, 1, 2, 3, 2),
                   (3, 2, 1, 3, 2, 3)]

    enhanced_chain = enhance_braid_move_chain(braid_chain,CartanType(["A",5])) 
    assert enhanced_chain[0][1] == None
    assert enhanced_chain[0][2] == None
    assert enhanced_chain[7][1] == (3, 6)
    assert enhanced_chain[7][2] == (-1, -1)


def diff_interval(list1,list2):
    """
    Return the smallest contiguous half-open interval [a,b) 
    that contains the indices where ``list1`` and ``list2`` differ. 
    Return ``None`` if the lists don't differ.

    INPUT:

    - ``list1``, ``list2`` -- two lists of the same length (no length checking)
    """
    min_diff = None
    max_diff = None
    for i,elem1 in enumerate(list1):
        if elem1 != list2[i]:
            min_diff = i
            break
    if min_diff is None:
        if len(list1) == len(list2):
            return None
        else:
            raise IndexError
    for j in xrange(len(list2), min_diff, -1):
        if list1[j-1] != list2[j-1]:
            max_diff = j
            break
    return (min_diff, max_diff)


def diff_interval_test():
    assert diff_interval([1,2,3,4],[1,2,3,4]) == None 
    assert diff_interval([1,2,3,4],[1,3,2,4]) == (1,3)
    assert diff_interval([1,2,4,4,5],[1,3,45,6,3]) == (1,5)
    try: 
        diff_interval([1],[2,3])
        raise AssertionError, "IndexError not raised"
    except IndexError:
        pass
    try: 
        diff_interval([],[2,3])
        raise AssertionError, "IndexError not raised"
    except IndexError:
        pass

