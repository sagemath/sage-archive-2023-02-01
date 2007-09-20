from sage.structure.sage_object import SageObject

class _ProofPref(SageObject):
    """
    An object that holds global proof preferences.  For now these are merely True/False flags for various parts of Sage that use probabilistic algorithms.
    A True flag means that the subsystem (such as linear algebra or number fields) should return results that are true unconditionally: the correctness should not depend on an algorithm with a nonzero probability of returning an incorrect answer or on the truth of any unproven conjectures.
    A False flag means that the subsystem can use faster methods to return answers that have a very small probability of being wrong.
    """
    def __init__(self, proof = True):
        self._require_proof = {}
        self._require_proof["arithmetic"] = proof
        self._require_proof["elliptic_curve"] = proof
        self._require_proof["linear_algebra"] = proof
        self._require_proof["number_field"] = proof
        self._require_proof["other"] = proof

    def arithmetic(self, t = None):
        """
        Controls the default proof strategy for integer arithmetic algorithms (such as primality testing).

        INPUT:
            t -- boolean or None

        OUTPUT:
            If t == True, requires integer arithmetic operations to (by default) return results that are true unconditionally: the correctness will not depend on an algorithm with a nonzero probability of returning an incorrect answer or on the truth of any unproven conjectures.
            If t == False, allows integer arithmetic operations to (by default) return results that may depend on unproven conjectures or on probabilistic algorithms.  Such algorithms often have a substantial speed improvement over those requiring proof.
            If t is None, returns the integer arithmetic proof status.

        EXAMPLES:
            sage: proof.arithmetic()
            True
            sage: proof.arithmetic(False)
            sage: proof.arithmetic()
            False
            sage: proof.arithmetic(True)
            sage: proof.arithmetic()
            True
        """
        if t is None:
            return self._require_proof["arithmetic"]
        self._require_proof["arithmetic"] = bool(t)

    def elliptic_curve(self, t = None):
        """
        Controls the default proof strategy for elliptic curve algorithms.

        INPUT:
            t -- boolean or None

        OUTPUT:
            If t == True, requires elliptic curve algorithms to (by default) return results that are true unconditionally: the correctness will not depend on an algorithm with a nonzero probability of returning an incorrect answer or on the truth of any unproven conjectures.
            If t == False, allows elliptic curve algorithms to (by default) return results that may depend on unproven conjectures or on probabilistic algorithms.  Such algorithms often have a substantial speed improvement over those requiring proof.
            If t is None, returns the current elliptic curve proof status.

        EXAMPLES:
            sage: proof.elliptic_curve()
            True
            sage: proof.elliptic_curve(False)
            sage: proof.elliptic_curve()
            False
            sage: proof.elliptic_curve(True)
            sage: proof.elliptic_curve()
            True
        """
        if t is None:
            return self._require_proof["elliptic_curve"]
        self._require_proof["elliptic_curve"] = bool(t)

    def linear_algebra(self, t = None):
        """
        Controls the default proof strategy for linear algebra algorithms.

        INPUT:
            t -- boolean or None

        OUTPUT:
            If t == True, requires linear algebra algorithms to (by default) return results that are true unconditionally: the correctness will not depend on an algorithm with a nonzero probability of returning an incorrect answer or on the truth of any unproven conjectures.
            If t == False, allows linear algebra algorithms to (by default) return results that may depend on unproven conjectures or on probabilistic algorithms.  Such algorithms often have a substantial speed improvement over those requiring proof.
            If t is None, returns the current linear algebra proof status.

        EXAMPLES:
            sage: proof.linear_algebra()
            True
            sage: proof.linear_algebra(False)
            sage: proof.linear_algebra()
            False
            sage: proof.linear_algebra(True)
            sage: proof.linear_algebra()
            True
        """
        if t is None:
            return self._require_proof["linear_algebra"]
        self._require_proof["linear_algebra"] = bool(t)

    def number_field(self, t = None):
        """
        Controls the default proof strategy for number field algorithms.

        INPUT:
            t -- boolean or None

        OUTPUT:
            If t == True, requires number field algorithms to (by default) return results that are true unconditionally: the correctness will not depend on an algorithm with a nonzero probability of returning an incorrect answer or on the truth of any unproven conjectures.
            If t == False, allows number field algorithms to (by default) return results that may depend on unproven conjectures or on probabilistic algorithms.  Such algorithms often have a substantial speed improvement over those requiring proof.
            If t is None, returns the current number field proof status.

        EXAMPLES:
            sage: proof.number_field()
            True
            sage: proof.number_field(False)
            sage: proof.number_field()
            False
            sage: proof.number_field(True)
            sage: proof.number_field()
            True
        """
        if t is None:
            return self._require_proof["number_field"]
        self._require_proof["number_field"] = bool(t)

_proof_prefs = _ProofPref(True) #Creates the global object that stores proof preferences.

def get_flag(t = None, subsystem = None):
    """
    Used for easily determining the correct proof flag to use.

    EXAMPLES:
        sage: from sage.structure.proof import get_flag
        sage: get_flag(None)
        True
        sage: get_flag(False)
        False
        sage: get_flag(True)
        True

    """
    if t is None:
        if subsystem in ["arithmetic", "elliptic_curve", "linear_algebra", "number_field"]:
            return _proof_prefs[subsystem]
        else:
            return _proof_prefs["other"]
    return t

