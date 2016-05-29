from operator import mul
from sage.matrix.constructor import matrix
from sage.functions.other import binomial
from sage.calculus.var import var
from sage.misc.functional import symbolic_sum
from sage.coding.linear_code import AbstractLinearCode, LinearCodeSyndromeDecoder
from sage.coding.encoder import Encoder
import sage.combinat as subsets
from sage.rings.finite_rings.finite_field_base import FiniteField

#to compute the sum of n chose i where i ranges from 0 to k
def binomialSum(n,k):
    s=1
    nCi=1
    for i in range(k):
        nCi=((n-i)*nCi)/(i+1)
        s=nCi+s
    return s

#to use the evaluation of polynomial at the points to obtain the polynomial 
def multivariatePolynomialInterpolation(evaluation, numberOfVariable, order, q, finiteField, _R):
    if numberOfVariable==0 or order==0:
        return evaluation[0]
    xcordinate=finiteField.list()
    nq=q**(numberOfVariable-1)
    d=min(order+1,q)
    evaluation2=[]
    uniPolyRing=PolynomialRing(finiteField,'x')
    for k in range(nq):
        points=[(xcordinate[i], evaluation[k+i*nq]) for i in range(q)]
        polyVector=uniPolyRing.lagrange_polynomial(points).coefficients(sparse=False)
        if len(polyVector)<d:
            #adding zeros to represet a (d-1) degree polynomial
            polyVector=polyVector+[0 for i in range(d-len(polyVector))]
        evaluation2.append(polyVector)
    poly=0
    z=1
    x=_R.gen(numberOfVariable-1)
    for k in range(d):#computing the polynomial
        poly=poly+z*multivariatePolynomialInterpolation([evaluation2[i][k] for i in range(nq)], numberOfVariable-1, order-k, q, finiteField, _R)
        z=z*x
    return poly

def ReedMullerCode(finiteField, order, numberOfVariable):
    if (isinstance(finiteField,FiniteField)):
            baseField=finiteField
            q=baseField.cardinality()
        elif (isinstance(finiteField, Integer)):
             baseField=GF(finiteField, 'x')
             q=finiteField
    if q == 2:
        return BinaryReedMullerCode(order, numberOfVariable)
    else:
        return QAryReedMullerCode(baseField, order, numberOfVariable)

class QAryReedMullerCode(AbstractLinearCode):
    _registered_encoders={}
    _registered_decoders={}

    def __init__(self, finiteField, order, numberOfVariable):
        #to handle the case when the base field is directly given and input sanitization
        if (isinstance(finiteField,FiniteField)):
            baseField=finiteField
            q=baseField.cardinality()
        elif (isinstance(finiteField, Integer)):
             baseField=GF(finiteField, 'x')
             q=finiteField
        else:
            raise ValueError("Incorrect data-type of input: You must either give the size of the finite field or the finite field itself")
        if not(isinstance(order,Integer)):
            raise ValueError("Incorrect data-type of input: The order of the code must be an integer")
        if not(isinstance(numberOfVariable, Integer)):
            raise ValueError("Incorrect data-type of input: The number of variables must be an integer")
        if (order>=q):
            raise ValueError("The order must be less than %s" % q)

        super(QAryReedMullerCode, self).__init__(baseField,q**numberOfVariable,"ReedMullerVectorEncoder","Syndrome")
        self.order=order
        self.numberOfVariable=numberOfVariable
        self.q=q
        self._dimension=binomial(numberOfVariable+order, order)

    def _repr_(self):
        return "%s-ary Reed Muller Code of order %s and number of variables %s" % (self.q, self.order, self.numberOfVariable)

    def _latex_(self):
        return "%s\textnormal{-ary Reed Muller Code of order} %s \textnormal{and number of variables} %s" % (self.q, self.order, self.numberOfVariable)

    def __eq__(self,other):
        return (isinstance(other, QAryReedMullerCode)) and self.q==other.q and self.order==other.order and self.numberOfVariable==other.numberOfVariable

class BinaryReedMullerCode(AbstractLinearCode):
    _registered_encoders={}
    _registered_decoders={}

    def __init__(self, order, numberOfVariable):
        #input sanitization
        if not(isinstance(order,Integer)):
            raise ValueError("Incorrect data-type of input: The order of the code must be an integer")
        if not(isinstance(numberOfVariable, Integer)):
            raise ValueError("Incorrect data-type of input: The number of variables must be an integer")
        if (numberOfVariable<order):
            raise ValueError("The order must be less than %s" % numberOfVariable)

        super(BinaryReedMullerCode, self).__init__(GF(2), 2**numberOfVariable,"ReedMullerVectorEncoder","Syndrome")
        self.order=order
        self.numberOfVariable=numberOfVariable
        self.q=2
        self._dimension=binomialSum(numberOfVariable,order)

    def _repr_(self):
        return "Binary Reed Muller Code of order %s and number of variables %s" % (self.order, self.numberOfVariable)

    def _latex_(self):
        return "\textnormal{Binary Reed Muller Code of order} %s \textnormal{and number of variables} %s" % (self.q, self.order, self.numberOfVariable)

    def __eq__(self,other):
        return (isinstance(other, BinaryReedMullerCode)) and self.order==other.order and self.numberOfVariable==other.numberOfVariable

class ReedMullerVectorEncoder(Encoder):
    def __init__(self, code):
        super(ReedMullerVectorEncoder, self).__init__(code)

    def _repr_(self):
        return "Evaluation polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        return "\textnormal{Evaluation polynomial-style encoder for }%s" % self.code()._latex_()

    def __eq__(self,other):
        return (isinstance(other, ReedMullerVectorEncoder)) and self.code==other.code

    def generator_matrix(self):
        C=self.code()
        baseField=C.base_field()
        order=C.order
        numberOfVariable=C.numberOfVariable
        q=C.q
        baseFieldTuple=Tuples(baseField.list(),numberOfVariable)
        exponents=Subsets(range(numberOfVariable)*(q-1), submultiset=True)[0:C.dimension()]
        return Matrix(baseField, [[reduce(mul, [x[i] for i in exponent],1) for x in baseFieldTuple] for exponent in exponents])

class ReedMullerPolynomialEncoder(Encoder):
    def __init__(self, code):
        super(ReedMullerPolynomialEncoder, self).__init__(code)
        self._R=PolynomialRing(self.code().base_field(),self.code().numberOfVariable,"x")

    def _repr_(self):
        return "Evaluation polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        return "\textnormal{Evaluation polynomial-style encoder for }%s" % self.code()._latex_()

    def __eq__(self,other):
        return (isinstance(other, ReedMullerPolynomialEncoder)) and self.code==other.code

    def encode(self, p):
        M = self.message_space()
        if p not in M:
            raise ValueError("The value to encode must be in %s" % M)
        C=self.code()
        if p.degree() > C.order:
            raise ValueError("The polynomial to encode must have degree at most %s" % C.order)
        baseFieldTuple = Tuples(C.base_field().list(), C.numberOfVariable)
        return vector(C.base_ring(), [p(x) for x in baseFieldTuple])
    
    def unencode_nocheck(self, c):
        return multivariatePolynomialInterpolation(c, self.code().numberOfVariable, self.code().order, self.code().q, self.code().base_field(), self._R)


    def message_space(self):
        return self._R

QAryReedMullerCode._registered_encoders["ReedMullerVectorEncoder"] = ReedMullerVectorEncoder
QAryReedMullerCode._registered_encoders["ReedMullerPolynomialEncoder"] = ReedMullerPolynomialEncoder

QAryReedMullerCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder

BinaryReedMullerCode._registered_encoders["ReedMullerVectorEncoder"] = ReedMullerVectorEncoder
BinaryReedMullerCode._registered_encoders["ReedMullerPolynomialEncoder"] = ReedMullerPolynomialEncoder

BinaryReedMullerCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
