#ifndef MINOR_H
#define MINOR_H

#include <assert.h>
#include <time.h>
#include <iostream>
#include <string>

using namespace std;

/*! \class MinorKey
    \brief Class MinorKey can be used for representing keys in a cache for
    sub-determinantes; see class Cache.

    As such, it is a realization of the template class KeyClass which is used
    in the declaration of class Cache. Following the documentation of class
    Cache, we need to implement at least the methods:<br>
    <c>bool MinorKey::operator< (const MinorKey& key),</c><br>
    <c>bool MinorKey::operator== (const MinorKey& key),</c><br>
    MinorKey uses two private arrays of ints \c _rowKey and \c _columnKey to
    encode rows and columns of a pre-defined matrix. Semantically, the row
    indices and column indices form the key for caching the value of the
    corresponding minor.<br>
    More concretely, let us assume that the pre-defined matrix has
    <em>32*R+r, r<32,</em> rows and <em>32*C+c, c<32,</em> columns. All row
    indices can then be captured using R+1 ints, since an int is a
    32-bit-number (regardless of the platform). The analog holds for the
    columns. Consequently, each instance of MinorKey encodes the sets of rows
    and columns which shall belong to the minor of interest (and which shall
    not).<br>
    Example: The \c _rowKey with \c _rowKey[1] = 0...011 and
    \c _rowKey[0] = 0...01101 encodes the rows with indices 33, 32, 3, 2,
    and 0.
    \author Frank Seelisch, http://www.mathematik.uni-kl.de/~seelisch
*/
class MinorKey
{
  private:
     /**
     * a pointer to an array[0..k-1] of ints, capturing k*32 bits for
     * determining which rows of a pre-defined matrix shall belong to the
     * minor of interest;
     * for i < j, _rowKey[i] holds lower bits than _rowKey[j]
     */
     unsigned int* _rowKey;

     /**
     * a pointer to an array[0..k-1] of ints, capturing k*32 bits for
     * determining which columns of a pre-defined matrix shall belong to the
     * minor of interest;
     * for i < j, _columnKey[i] holds lower bits than _columnKey[j]
     */
     unsigned int* _columnKey;

     /**
     * the number of ints (i.e. 32-bit-numbers) we need to encode the set of
     * rows;
     * If the higest row index is 70, we need 3 blocks of 32 bits to also
     * encode the 70th bit.
     */
     int _numberOfRowBlocks;

     /**
     * the number of ints (i.e. 32-bit-numbers) we need to encode the set of
     * columns;
     * If the higest column index is 70, we need 3 blocks of 32 bits to also
     * encode the 70th bit.
     */
     int _numberOfColumnBlocks;

     /**
     * Inlined accessor of blockIndex-th element of _rowKey.
     * @param blockIndex the index of the int to be retrieved
     * @return an entry of _rowKey
     */
     unsigned int getRowKey (const int blockIndex) const;

     /**
     * Accessor of blockIndex-th element of _columnKey.
     * @param blockIndex the index of the int to be retrieved
     * @return an entry of _columnKey
     */
     unsigned int getColumnKey (const int blockIndex) const;

     /**
     * A method for setting the blockIndex-th element of _rowKey.
     * @param blockIndex the index of the int to be retrieved
     * @param rowKey the row key to be set
     */
     void setRowKey (const int blockIndex, const unsigned int rowKey);

     /**
     * A method for setting the blockIndex-th element of _columnKey.
     * @param blockIndex the index of the int to be retrieved
     * @param columnKey the column key to be set
     */
     void setColumnKey (const int blockIndex, const unsigned int columnKey);

     /**
     * Accessor of _numberOfRowBlocks.
     * @return the number of 32-bit-blocks needed to encode all rows of the
     *         minor as a sequence of bits
     */
     int getNumberOfRowBlocks () const;

     /**
     * Accessor of _numberOfColumnBlocks.
     * @return the number of 32-bit-blocks needed to encode all columns of
     *         the minor as a sequence of bits
     */
     int getNumberOfColumnBlocks () const;

     /**
     * A method for deleting all entries of _rowKey and _columnKey.
     */
     void reset();

     /**
     * A method for counting the number of set bits.
     * For a == 1, the number of set bits in _rowKey will be counted;
     * for a == 2 in _columnKey.
     * This method will only be called in the debug version.
     * @param a for controlling whether to count in _rowKey or _columnKey
     * @return the number of set bits either in _rowKey or _columnKey
     */
     int getSetBits (const int a) const;

     /**
     * For letting MinorProcessor see the private methods of this class.
     */
     friend class MinorProcessor;
  public:
     /**
     * A constructor for class MinorKey.
     * The ints given in the array rowKey encode all rows which shall belong
     * to the minor. Each array entry encodes 32 rows, e.g. the i-th array
     * entry 0...01101 encodes the rows with absolute matrix row indices
     * 3+i*32, 2+i*32, and 0+i*32. Analog for columns.
     * @param lengthOfRowArray the length of the array rowKey
     * @param rowKey a pointer to an array of ints encoding the set of rows of
     *        the minor
     * @param lengthOfColumnArray the length of the array columnKey
     * @param columnKey a pointer to an array of ints encoding the set of
     +        columns of the minor
     */
     MinorKey (const int lengthOfRowArray = 0,
               const unsigned int* const rowKey = 0,
               const int lengthOfColumnArray = 0,
               const unsigned int* const columnKey = 0);

     /**
     * A setter method for class MinorKey.
     * Just like the constructor of this class, this method will set all
     * private fields according to the given parameters. Note that this method
     * will change the given instance of MinorKey.
     * @param lengthOfRowArray the length of the array rowKey
     * @param rowKey a pointer to an array of ints encoding the set of rows of
     * the minor
     * @param lengthOfColumnArray the length of the array columnKey
     * @param columnKey a pointer to an array of ints encoding the set of
     *        columns of the minor
     * @see MinorKey::MinorKey (const int lengthOfRowArray, const int* rowKey,
                                const int lengthOfColumnArray,
                                const int* columnKey)
     */
     void set(const int lengthOfRowArray, const unsigned int* rowKey,
              const int lengthOfColumnArray, const unsigned int* columnKey);

     /**
     * A copy constructor.
     * This method overrides the shallow copy constructor by a self-written
     * deep copy version.
     * @param mk the MinorKey to be deep copied
     */
     MinorKey (const MinorKey& mk);

     /**
     * A destructor for deleting an instance.
     */
     ~MinorKey ();

     /**
     * just to make the compiler happy
     */
     MinorKey& operator=(const MinorKey&);

     /**
     * just to make the compiler happy
     */
     bool operator==(const MinorKey&) const;

     /**
     * just to make the compiler happy
     */
     bool operator<(const MinorKey&) const;

     /**
     * A method for retrieving the (0-based) index of the i-th row in the set
     * of rows encoded in \a this.
     * Lower values for \c i result in lower absolute row indices.
     * \par Example:
     * Applied to the row pattern 10010001101 and i = 3, we get the 0-based
     * index of the 3-rd set bit counted from the right, i.e. 7.
     * \par Assertion
     * The method assumes that there are at least \c i rows encoded in the
     * given MinorKey.
     * @param i the relative index of the row, as encoded in \a this
     * @return (0-based) absolute row index of the i-th row in \a this
     */
     int getAbsoluteRowIndex (const int i) const;

     /**
     * A method for retrieving the (0-based) index of the i-th column in the
     * set of columns encoded in \a this.
     * Lower values for \c i result in lower absolute column indices.
     * \par Example:
     * Applied to the column pattern 10010001101 and i = 3, we get the 0-based
     * index of the 3-rd set bit counted from the right, i.e. 7.
     * \par Assertion
     * The method assumes that there are at least \c i columns encoded in the
     * given MinorKey.
     * @param i the relative index of the column, as encoded in \a this
     * @return (0-based) absolute column index of the i-th row in \a this
     */
     int getAbsoluteColumnIndex (const int i) const;

     /**
     * A method for retrieving the (0-based) relative index of the i-th row
     * in \a this MinorKey.
     * Lower values for \c i result in lower relative row indices. Note that
     * the absolute index \c i is 0-based, too.
     * \par Example:
     * Applied to the row pattern 10010001101 and i = 7, we get the relative
     * 0-based position of the bit representing the absolute index 7, i.e. 3.
     * \par Assertion
     * The method assumes that the bit which corresponds to the absolute index
     * i is actually set.
     * @param i the absolute 0-based index of a row encoded in \a this
     * @return (0-based) relative row index corresponding to \c i
     */
     int getRelativeRowIndex (const int i) const;

     /**
     * A method for retrieving the (0-based) relative index of the i-th column
     * in \a this MinorKey.
     * Lower values for \c i result in lower relative column indices. Note that
     * the absolute index \c i is 0-based, too.
     * \par Example:
     * Applied to the column pattern 10010001101 and i = 7, we get the
     * relative 0-based position of the bit representing the absolute index 7,
     * i.e. 3.
     * \par Assertion
     * The method assumes that the bit which corresponds to the absolute index
     * i is actually set.
     * @param i the absolute 0-based index of a column encoded in \a this
     * @return (0-based) relative column index corresponding to \c i
     */
     int getRelativeColumnIndex (const int i) const;

     /**
     * A method for retrieving the 0-based indices of all rows encoded in \a
     * this MinorKey.
     * The user of this method needs to know the number of rows in \a this,
     * in order to know which indices in \c target[k] will be valid.
     * \par Example:
     * The bit pattern <c>0...01101</c> will give rise to the settings
     * <c>target[0] = 0, target[1] = 2, target[2] = 3</c>, and the user needs
     * to know in advance that there are three rows in \a this MinorKey.
     * \par Assertion
     * The method assumes that target has enough positions for all rows
     * encoded in \a this MinorKey.
     * @param target a pointer to some array of ints that is to be filled with
     *        the requested indices
     */
     void getAbsoluteRowIndices(int* const target) const;

     /**
     * A method for retrieving the 0-based indices of all columns encoded in
     * \a this MinorKey.
     * The user of this method needs to know the number of columns in \a this,
     * in order to know which indices in \c target[k] will be valid.
     * \par Example:
     * The bit pattern <c>0...01101</c> will give rise to the settings
     * <c>target[0] = 0, target[1] = 2, target[2] = 3</c>, and the user needs
     * to know in advance that there are three columns in \a this MinorKey.
     * \par Assertion
     * The method assumes that target has enough positions for all columns
     * encoded in \a this MinorKey.
     * @param target a pointer to some array of ints that is to be filled
     *        with the requested indices
     */
     void getAbsoluteColumnIndices(int* const target) const;

     /**
     * A method for retrieving a sub-MinorKey resulting from omitting one row
     * and one column of \a this MinorKey.
     * \par Assertion
     * The method assumes that the row with absolute index
     * \c absoluteEraseRowIndex (counted from lower bits to higher bits) and
     * the column with absolute index \c absoluteEraseColumnIndex are actually
     * set in \c mk.
     * @param absoluteEraseRowIndex the 0-based absolute index of a row in
     *        \a mk
     * @param absoluteEraseColumnIndex the 0-based absolute index of a column
     *        in \a mk
     * @return the MinorKey when omitting the specified row and column
     */
     MinorKey getSubMinorKey (const int absoluteEraseRowIndex,
                              const int absoluteEraseColumnIndex) const;

     /**
     * A comparator for two instances of MinorKey.
     * The ordering induced by this implementation determines the ordering of
     * all (key --> value) pairs in a cache that uses MinorKey as KeyClass.
     * @param mk a second MinorKey to be compared with \a this instance
     * @return -1 iff \a this instance is smaller than \a mk; 0 for equality;
     *         +1 otherwise
     * @see MinorKey::operator== (const MinorKey&) const
     */
     int compare (const MinorKey& mk) const;

     /**
     * This method redefines the set of rows represented by \a this MinorKey.
     * After the method, the defined set of rows coincides with the lowest
     * \c k rows of \c mk. (Here, lowest means w.r.t. indices.)<br>
     * Note that the method modifies the given instance of MinorKey.
     * \par Assertion
     * It is assumed that \c mk represents at least \c k rows.
     * @param k the number of rows
     * @param mk the MinorKey from which to choose the lowest \c k rows
     * @see MinorKey::selectNextRows (const int k, const MinorKey& mk)
     */
     void selectFirstRows (const int k, const MinorKey& mk);

     /**
     * This method redefines the set of rows represented by \a this MinorKey.
     * Both the old and the new set of \c k rows are subsets of the rows
     * represented by \c mk. After the method, the defined set of rows is
     * the next sensible choice of \c k rows of \c mk. (Here, next means
     * the next w.r.t. the increasing index ordering on multi-indices of
     * natural numbers.)<br>
     * Note that the method modifies the given instance of MinorKey.
     * \par Assertion
     * It is assumed that \c mk represents at least \c k rows. Furthermore,
     * the method assumes that the old set of rows represented by \a this
     * is also a subset of the rows given by \c mk.
     * @param k the number of rows
     * @param mk the MinorKey from which to choose the lowest \c k rows
     * @return true iff there is a next choice of \c k rows
     * @see MinorKey::selectFirstRows (const int k, const MinorKey& mk)
     */
     bool selectNextRows (const int k, const MinorKey& mk);

     /**
     * This method redefines the set of columns represented by \a this
     * MinorKey.
     * After the method, the defined set of columns coincides with the lowest
     * \c k columns of \c mk. (Here, lowest means w.r.t. indices.)<br>
     * Note that the method modifies the given instance of MinorKey.
     * \par Assertion
     * It is assumed that \c mk represents at least \c k columns.
     * @param k the number of columns
     * @param mk the MinorKey from which to choose the lowest \c k columns
     * @see MinorKey::selectNextColumns (const int k, const MinorKey& mk)
     */
     void selectFirstColumns (const int k, const MinorKey& mk);

     /**
     * This method redefines the set of columns represented by \a this
     * MinorKey.
     * Both the old and the new set of \c k columns are subsets of the columns
     * represented by \c mk. After the method, the defined set of columns is
     * the next sensible choice of \c k columns of \c mk. (Here, next means
     * the next w.r.t. the increasing index ordering on multi-indices of
     * natural numbers.)<br>
     * Note that the method modifies the given instance of MinorKey.
     * \par Assertion
     * It is assumed that \c mk represents at least \c k columns. Furthermore,
     * the method assumes that the old set of columns represented by \a this
     * is also a subset of the columns given by \c mk.
     * @param k the number of columns
     * @param mk the MinorKey from which to choose the lowest \c k columns
     * @return true iff there is a next choice of \c k columns
     * @see MinorKey::selectFirstColumns (const int k, const MinorKey& mk)
     */
     bool selectNextColumns (const int k, const MinorKey& mk);

     /**
     * A method for providing a printable version of the represented MinorKey.
     * @return a printable version of the given instance as instance of class
     * string
     */
     string toString () const;

     /**
     * A method for printing a string representation of the given MinorKey to
     * std::cout.
     */
     void print () const;
};

class MinorValue
{
  protected:
    /**
    * -1 iff cache is not used, otherwise the number of retrievals so far of
    * the current minor
    */
    int _retrievals;

    /**
    * -1 iff cache is not used, otherwise the maximum number of potential
    * retrievals of this minor (e.g. when the minor would be kept in cache
    * forever)
    */
    int _potentialRetrievals;

    /**
    * a store for the actual number of multiplications to compute the current
    * minor
    */
    int _multiplications;

    /**
    * a store for the actual number of additions to compute the current minor
    */
    int _additions;

    /**
    * a store for the accumulated number of multiplications to compute the
    * current minor;
    * This also includes all multiplications nested in sub-minors which may be
    * retrieved from a cache. (Thus, these nested operations do not need to be
    * performed again.)
    */
    int _accumulatedMult;

    /**
    * a store for the accumulated number of additions to compute the current
    * minor;
    * This also includes all additions nested in sub-minors which may be
    * retrieved from a cache. (Thus, these nested operations do not need to be
    * performed again.)
    */
    int _accumulatedSum;

    /**
    * A method for obtaining a rank measure for the given MinorValue.<br>
    * Rank measures are used to compare any two instances of MinorValue. The
    * induced ordering
    * on MinorValues has an impact on the caching behaviour in a given cache:
    * Greater MinorValues will be cached longer than lower ones.<br>
    * More explicitely, this means: Make the return value of this method
    * greater, and the given MinorValue will be cached longer when caching
    * strategy 1 is deployed.<br>
    * Rank measure 1 is equal to the number of actually performed
    * multiplications to compute \a mv.
    * @return an integer rank measure of \c this
    * @see MinorValue::operator< (const MinorValue& mv)
    */
    int rankMeasure1 () const;

    /**
    * A method for obtaining a rank measure for the given MinorValue.<br>
    * Rank measures are used to compare any two instances of MinorValue. The
    * induced ordering on MinorValues has an impact on the caching behaviour
    * in a given cache: Greater MinorValues will be cached longer than lower
    * ones.<br>
    * More explicitely, this means: Make the return value of this method
    * greater, and the given MinorValue will be cached longer when caching
    * strategy 1 is deployed.<br>
    * Rank measure 2 is equal to the number of accumulated multiplications to
    * compute the given MinorValue. This also includes all nested
    * multiplications which were performed to compute all sub-minors which
    * could be reused from cache.
    * @return an integer rank measure of \c this
    * @see MinorValue::operator< (const MinorValue& mv)
    */
    int rankMeasure2 () const;

    /**
    * A method for obtaining a rank measure for the given MinorValue.<br>
    * Rank measures are used to compare any two instances of MinorValue. The
    * induced ordering on MinorValues has an impact on the caching behaviour
    * in a given cache: Greater MinorValues will be cached longer than lower
    * ones.<br>
    * More explicitely, this means: Make the return value of this method
    * greater, and the given MinorValue will be cached longer when caching
    * strategy 1 is deployed.<br>
    * Rank measure 3 is equal to the number of actually performed
    * multiplications, weighted with the ratio of not yet performed retrievals
    * over the maximum number of retrievals.
    * @return an integer rank measure of \c this
    * @see MinorValue::operator< (const MinorValue& mv)
    */
    int rankMeasure3 () const;

    /**
    * A method for obtaining a rank measure for the given MinorValue.<br>
    * Rank measures are used to compare any two instances of MinorValue. The
    * induced ordering on MinorValues has an impact on the caching behaviour
    * in a given cache: Greater MinorValues will be cached longer than lower
    * ones.<br>
    * More explicitely, this means: Make the return value of this method
    * greater, and the given MinorValue will be cached longer when caching
    * strategy 1 is deployed.<br>
    * Rank measure 4 is equal to the number of actually performed
    * multiplications, multiplied with the number of not yet performed
    * retrievals.
    * @return an integer rank measure of \c this
    * @see MinorValue::operator< (const MinorValue& mv)
    */
    int rankMeasure4 () const;

    /**
    * A method for obtaining a rank measure for the given MinorValue.<br>
    * Rank measures are used to compare any two instances of MinorValue. The
    * induced ordering on MinorValues has an impact on the caching behaviour
    * in a given cache: Greater MinorValues will be cached longer than lower
    * ones.<br>
    * More explicitely, this means: Make the return value of this method
    * greater, and the given MinorValue will be cached longer when caching
    * strategy 1 is deployed.<br>
    * Rank measure 5 is equal to the number of not yet performed retrievals.
    * This strategy tends to cache MinorValues longer which have a high
    * maximum number of potential retrievals.
    * @return an integer rank measure of \c this
    * @see MinorValue::operator< (const MinorValue& mv)
    */
    int rankMeasure5 () const;

    /**
    * private store for the current value ranking strategy;
    * This member can be set using MinorValue::SetRankingStrategy (const int).
    */
    static int _RankingStrategy;

    /**
    * Accessor for the static private field _RankingStrategy.
    */
    static int GetRankingStrategy();
  public:
    /**
    * just to make the compiler happy
    */
    bool operator== (const MinorValue& mv) const;

    /**
    * just to make the compiler happy
    */
    bool operator< (const MinorValue& mv) const;

    /**
    * A method for retrieving the weight of a given MinorValue.
    * The implementation of Cache uses this function to determine the total
    * weight of an entire cache. As the user can instantiate Cache by
    * determining its maximum total weight
    * (see Cache::Cache(const int, const int)),
    * the definition of weight of a MinorValue
    * may have an impact on the behaviour of the cache.
    * @return the weight of a given instance of MinorValue
    * @see Cache::getWeight () const
    */
    virtual int getWeight () const;

    /**
    * A method for accessing the number of retrievals of this minor. Multiple
    * retrievals will occur when computing large minors by means of cached
    * sub-minors. (Then, the latter ones may be retrieved multiple times.)
    * @return the number of retrievals of this minor
    * @see MinorValue::getPotentialRetrievals () const
    */
    int getRetrievals () const;

    /**
    * A method for accessing the maximum number of potential retrievals of
    * this minor. Multiple retrievals will occur when computing large minors
    * by means of cached sub-minors. (Then, the latter ones may be retrieved
    * multiple times.)
    * @return the maximum number of potential retrievals of this minor
    * @see MinorValue::getRetrievals () const
    */
    int getPotentialRetrievals () const;

    /**
    * A method for accessing the multiplications performed while computing
    * this minor.
    * Due to multiplication with zero entries of the underlying matrix, some
    * sub-minors may be irrelevant. In this case, the multiplications needed
    * to compute these sub-minors will not be counted (, as they need not be
    * performed).
    * Moreover, multiplications that were needed to compute cached sub-minors
    * will not be counted either, as the value of those sub-minors can be
    * directly retrieved from the cache.
    * @return the number of multiplications performed
    * @see MinorValue::getAccumulatedMultiplications () const
    */
    int getMultiplications () const;

    /**
    * A method for accessing the multiplications performed while computing
    * this minor, including all nested multiplications.
    * Contrary to MinorValue::getMultiplications () const, this method will
    * also count multiplications needed to compute all cached sub-minors
    * (, although they need not be performed again in order to compute the
    * given instance of MinorValue).
    * @return the number of multiplications performed, including nested
    *         multiplications
    * @see MinorValue::getMultiplications () const
    */
    int getAccumulatedMultiplications () const;

    /**
    * A method for accessing the additions performed while computing this
    * minor.
    * Additions that were needed to compute cached sub-minors will not be
    * counted, as the value of those sub-minors can be directly retrieved
    * from the cache.
    * @return the number of additions performed
    * @see MinorValue::getAccumulatedAdditions () const
    */
    int getAdditions () const;

    /**
    * A method for accessing the additions performed while computing this
    * minor, including all nested additions.
    * Contrary to MinorValue::getAdditions () const, this method will also
    * count additions needed to compute all cached sub-minors (, although
    * they need not be performed again in order to compute the given instance
    * of MinorValue).
    * @return the number of additions performed, including nested additions
    * @see MinorValue::getAdditions () const
    */
    int getAccumulatedAdditions () const;

    /**
    * A method for incrementing the number of performed retrievals of \a this
    * instance of MinorValue.<br>
    * Note that, when calling MinorValue::incrementRetrievals () for some
    * instance \a mv of MinorValue which has been cached in a Cache under
    * MinorKey \a mk, the user should be careful: After incrementing the
    * number of retrievals for \a mv, the user should always put the value
    * again into cache, i.e. should perform
    * Cache::put (const KeyClass&, const ValueClass&)
    * with \a mk and the modified \a mv as arguments. This is due to the fact
    * that changing the number of performed retrievals of a MinorValue may
    * have an impact on its ranking in Cache. Only by calling
    * Cache::put (const KeyClass&, const ValueClass&) can the user ensure
    * that the pair (\a mk --> \a mv) will be correctly re-positioned within
    * the Cache.
    */
    void incrementRetrievals ();

    /**
    * A method for obtaining a rank measure for theiven MinorValue.<br>
    * Rank measures are used to compare any two instances of MinorValue. The
    * induced ordering on MinorValues has an impact on the caching behaviour
    * of the underlying cache: Greater MinorValues will be cached longer than
    * lower ones.<br>
    * More explicitely, this means: Make the return value of this method
    * greater, and the given MinorValue will be cached longer.<br>
    * Internally, this method will call one of several implementations,
    * depending on the pre-defined caching strategy; see
    * MinorProcessor::SetCacheStrategy (const int).
    * @return an integer rank measure of \c this
    * @see MinorValue::operator< (const MinorValue& mv)
    * @see MinorProcessor::SetCacheStrategy (const int)
    */
    int getUtility () const;

    /**
    * A method for determining the value ranking strategy.<br>
    * This setting has a direct effect on how long the given MinorValue
    * will be cached in any cache that uses MinorValue to represent its
    * cached values.
    * @param rankingStrategy an int, so far one of 1, 2, ..., 5
    */
    static void SetRankingStrategy (const int rankingStrategy);

    /**
    * A method for providing a printable version of the represented MinorValue.
    * @return a printable version of the given instance as instance of class
    *         string
    */
    virtual string toString () const;

    /**
    * A method for printing a string representation of the given MinorValue
    * to std::cout.
    */
    void print () const;
};

/*! \class IntMinorValue
    \brief Class IntMinorValue is derived from MinorValue and can be used for
    representing values in a cache for sub-determinantes; see class Cache.

    As such, it is a realization of the template class ValueClass which is
    used in the declaration of class Cache. Following the documentation of
    class Cache, we need to implement at least the methods:<br>
    <c>bool IntMinorValue::operator< (const IntMinorValue& key),</c><br>
    <c>bool IntMinorValue::operator== (const IntMinorValue& key),</c><br>
    <c>int IntMinorValue::getWeight ().</c><br><br>
    The main purpose of IntMinorValue is to represent values of
    sub-determinantes of a pre-defined matrix. Class MinorKey is used to
    determine which rows and columns of this pre-defined matrix ought to
    belong to the respective sub-determinante of interest. So far,
    IntMinorValue is just an example implementation which assumes matrices
    with int entries, such that the result of any minor is again an int.
    \author Frank Seelisch, http://www.mathematik.uni-kl.de/~seelisch
*/
class IntMinorValue : public MinorValue
{
  private:
    /**
    * a store for the actual value of the minor
    */
    int _result;
  public:
    /**
    * A constructor for class MinorValue.
    * @param result the actual value of the represented minor
    * @param multiplications number of multiplications to compute \a this
             minor
    * @param additions number of additions to compute \a this minor
    * @param accumulatedMultiplications number of multiplications to compute
             \a this minor, including nested operations
    * @param accumulatedAdditions number of additions to compute \a this minor,
             including nested operations
    * @param retrievals number of times this minor has been retrieved from
             cache
    * @param potentialRetrievals maximum number of times this minor may be
             retrieved from cache
    */
    IntMinorValue (const int result, const int multiplications,
                   const int additions,
                   const int accumulatedMultiplications,
                   const int accumulatedAdditions, const int retrievals,
                   const int potentialRetrievals);

    /**
    * Copy constructor
    */
    IntMinorValue (const IntMinorValue& mv);

    /**
    * just to make the compiler happy
    */
    IntMinorValue ();

    /**
    * Destructor
    */
    virtual ~IntMinorValue ();

    /**
    * Accessor for the private field _result.
    * @result the result encoded in this class instance
    */
    int getResult() const;

    /**
    * Accessor for the current weight of this class instance.
    * @result the current weight of this class instance
    */
    int getWeight () const;

    /**
    * A method for providing a printable version of the represented MinorValue.
    * @return a printable version of the given instance as instance of class
    * string
    */
    string toString () const;
};

/*! \class PolyMinorValue
    \brief Class PolyMinorValue is derived from MinorValue and can be used for
    representing values in a cache for sub-determinantes; see class Cache.

    As such, it is a realization of the template class ValueClass which is
    used in the declaration of class Cache. Following the documentation of
    class Cache, we need to implement at least the methods:<br>
    <c>bool IntMinorValue::operator< (const IntMinorValue& key),</c><br>
    <c>bool IntMinorValue::operator== (const IntMinorValue& key),</c><br>
    <c>int IntMinorValue::getWeight ().</c><br><br>
    The main purpose of PolyMinorValue is to represent values of
    sub-determinantes of a pre-defined matrix. Class MinorKey is used to
    determine which rows and columns of this pre-defined matrix ought to
    belong to the respective sub-determinante of interest. PolyMinorValue is
    a special implementation which assumes matrices with polynomial entries,
    such that the result of any minor is again a polynomial.
    \author Frank Seelisch, http://www.mathematik.uni-kl.de/~seelisch
*/
class PolyMinorValue : public MinorValue
{
  private:
    /**
    * a store for the actual value of the minor
    */
    poly _result;
  public:
    /**
    * A constructor for class MinorValue.
    * @param result the actual value of the represented minor
    * @param multiplications number of multiplications to compute \a this
             minor
    * @param additions number of additions to compute \a this minor
    * @param accumulatedMultiplications number of multiplications to compute
             \a this minor, including nested operations
    * @param accumulatedAdditions number of additions to compute \a this
             minor, including nested operations
    * @param retrievals number of times this minor has been retrieved from
             cache
    * @param potentialRetrievals maximum number of times this minor may be
             retrieved from cache
    */
    PolyMinorValue (const poly result, const int multiplications,
                    const int additions,
                    const int accumulatedMultiplications,
                    const int accumulatedAdditions, const int retrievals,
                    const int potentialRetrievals);

    /**
    * Copy constructor for creating a deep copy.
    */
    PolyMinorValue (const PolyMinorValue& mv);

    /**
    * Assignment operator which creates a deep copy.
    */
    void operator= (const PolyMinorValue& mv);

    /**
    * just to make the compiler happy
    */
    PolyMinorValue ();

    /**
    * Destructor
    */
    virtual ~PolyMinorValue ();

    /**
    * Accessor for the private field _result.
    * @result the result encoded in this class instance
    */
    poly getResult() const;

    /**
    * Accessor for the current weight of this class instance.
    * @result the current weight of this class instance
    */
    int getWeight () const;

    /**
    * A method for providing a printable version of the represented MinorValue.
    * @return a printable version of the given instance as instance of class
    * string
    */
    string toString () const;
};

#endif
/* MINOR_H */
