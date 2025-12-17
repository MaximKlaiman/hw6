#ifndef HT_H
#define HT_H

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <functional>

typedef size_t HASH_INDEX_T;


// =======================
// Base Prober
// =======================
template <typename KeyType>
struct Prober {
    HASH_INDEX_T start_;
    HASH_INDEX_T m_;
    size_t numProbes_;
    static const HASH_INDEX_T npos = (HASH_INDEX_T)-1;

    void init(HASH_INDEX_T start, HASH_INDEX_T m, const KeyType& key)
    {
        (void) key;
        start_ = start;
        m_ = m;
        numProbes_ = 0;
    }

    HASH_INDEX_T next() {
        throw std::logic_error("Not implemented");
    }
};


// =======================
// Linear Prober
// =======================
template <typename KeyType>
struct LinearProber : public Prober<KeyType> {

    HASH_INDEX_T next()
    {
        if(this->numProbes_ >= this->m_) {
            return this->npos;
        }
        HASH_INDEX_T loc =
            (this->start_ + this->numProbes_) % this->m_;
        this->numProbes_++;
        return loc;
    }
};


// =======================
// Double Hash Prober
// =======================
template <typename KeyType, typename Hash2>
struct DoubleHashProber : public Prober<KeyType>
{
    Hash2 h2_;
    HASH_INDEX_T dhstep_;

    static const HASH_INDEX_T DOUBLE_HASH_MOD_VALUES[];
    static const int DOUBLE_HASH_MOD_SIZE;

private:
    HASH_INDEX_T findModulusToUseFromTableSize(HASH_INDEX_T currTableSize)
    {
        HASH_INDEX_T mod = DOUBLE_HASH_MOD_VALUES[0];
        for(int i = 0; i < DOUBLE_HASH_MOD_SIZE &&
            DOUBLE_HASH_MOD_VALUES[i] < currTableSize; i++)
        {
            mod = DOUBLE_HASH_MOD_VALUES[i];
        }
        return mod;
    }

public:
    DoubleHashProber(const Hash2& h2 = Hash2()) : h2_(h2) {}

    void init(HASH_INDEX_T start, HASH_INDEX_T m, const KeyType& key)
    {
        Prober<KeyType>::init(start, m, key);
        HASH_INDEX_T mod = findModulusToUseFromTableSize(m);
        dhstep_ = mod - (h2_(key) % mod);
    }

    HASH_INDEX_T next()
    {
        if(this->numProbes_ >= this->m_) {
            return this->npos;
        }
        HASH_INDEX_T loc =
            (this->start_ + this->numProbes_ * dhstep_) % this->m_;
        this->numProbes_++;
        return loc;
    }
};


// Static init for double hash
template <typename KeyType, typename Hash2>
const HASH_INDEX_T DoubleHashProber<KeyType, Hash2>::DOUBLE_HASH_MOD_VALUES[] =
{
    7, 19, 43, 89, 193, 389, 787, 1583, 3191, 6397, 12841, 25703, 51431, 102871,
    205721, 411503, 823051, 1646221, 3292463, 6584957, 13169963, 26339921,
    52679927, 105359939, 210719881, 421439749, 842879563, 1685759113
};

template <typename KeyType, typename Hash2>
const int DoubleHashProber<KeyType, Hash2>::DOUBLE_HASH_MOD_SIZE =
    sizeof(DoubleHashProber<KeyType, Hash2>::DOUBLE_HASH_MOD_VALUES) /
    sizeof(HASH_INDEX_T);


// =======================
// Hash Table
// =======================
template<
    typename K,
    typename V,
    typename Prober = LinearProber<K>,
    typename Hash = std::hash<K>,
    typename KEqual = std::equal_to<K> >
class HashTable
{
public:
    typedef K KeyType;
    typedef V ValueType;
    typedef std::pair<KeyType, ValueType> ItemType;
    typedef Hash Hasher;

    struct HashItem {
        ItemType item;
        bool deleted;
        HashItem(const ItemType& it) : item(it), deleted(false) {}
    };

    HashTable(double resizeAlpha = 0.4,
              const Prober& prober = Prober(),
              const Hasher& hash = Hasher(),
              const KEqual& kequal = KEqual());

    ~HashTable();

    bool empty() const;
    size_t size() const;

    void insert(const ItemType& p);
    void remove(const KeyType& key);

    ItemType const * find(const KeyType& key) const;
    ItemType * find(const KeyType& key);

    const ValueType& at(const KeyType& key) const;
    ValueType& at(const KeyType& key);
    const ValueType& operator[](const KeyType& key) const;
    ValueType& operator[](const KeyType& key);

    void reportAll(std::ostream& out) const;
    void clearTotalProbes() { totalProbes_ = 0; }
    size_t totalProbes() const { return totalProbes_; }

private:
    HashItem* internalFind(const KeyType& key) const;
    HASH_INDEX_T probe(const KeyType& key) const;
    void resize();

    static const HASH_INDEX_T npos = Prober::npos;
    static const HASH_INDEX_T CAPACITIES[];

    std::vector<HashItem*> table_;
    Hasher hash_;
    KEqual kequal_;
    mutable Prober prober_;
    mutable size_t totalProbes_;

    HASH_INDEX_T mIndex_;
    double resizeAlpha_;
    size_t currSize_;

    
    size_t occupied_;
};


// Capacities
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
const HASH_INDEX_T HashTable<K,V,Prober,Hash,KEqual>::CAPACITIES[] =
{
    11, 23, 47, 97, 197, 397, 797, 1597, 3203, 6421, 12853, 25717, 51437,
    102877, 205759, 411527, 823117, 1646237, 3292489, 6584983, 13169977,
    26339969, 52679969, 105359969, 210719881, 421439783, 842879579, 1685759167
};


// =======================
// HashTable Implementation
// =======================
template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
HashTable<K,V,Prober,Hash,KEqual>::HashTable(
    double resizeAlpha, const Prober& prober,
    const Hasher& hash, const KEqual& kequal)
    : hash_(hash), kequal_(kequal), prober_(prober)
{
    mIndex_ = 0;
    table_.resize(CAPACITIES[mIndex_], nullptr);
    resizeAlpha_ = resizeAlpha;
    currSize_ = 0;
    occupied_ = 0;
    totalProbes_ = 0;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
HashTable<K,V,Prober,Hash,KEqual>::~HashTable()
{
    for(size_t i = 0; i < table_.size(); i++) {
        if(table_[i] != nullptr) {
            delete table_[i];
        }
    }
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
bool HashTable<K,V,Prober,Hash,KEqual>::empty() const
{
    return currSize_ == 0;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
size_t HashTable<K,V,Prober,Hash,KEqual>::size() const
{
    return currSize_;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::insert(const ItemType& p)
{
    
    double lf = (double)occupied_ / (double)CAPACITIES[mIndex_];
    if(lf >= resizeAlpha_) {
        resize();
    }

    HASH_INDEX_T loc = probe(p.first);
    if(loc == npos) {
        throw std::logic_error("No space");
    }

    if(table_[loc] == nullptr) {
        table_[loc] = new HashItem(p);
        currSize_++;
        occupied_++; 
    }
    else {
        table_[loc]->item.second = p.second;
        if(table_[loc]->deleted) {
            table_[loc]->deleted = false;
            currSize_++;
            
        }
    }
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::remove(const KeyType& key)
{
    HASH_INDEX_T loc = probe(key);
    if(loc == npos) return;

    if(table_[loc] != nullptr && !table_[loc]->deleted) {
        table_[loc]->deleted = true;
        currSize_--;
        
    }
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::ItemType const *
HashTable<K,V,Prober,Hash,KEqual>::find(const KeyType& key) const
{
    HASH_INDEX_T h = probe(key);
    if(h == npos || table_[h] == nullptr) return nullptr;
    return &table_[h]->item;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::ItemType *
HashTable<K,V,Prober,Hash,KEqual>::find(const KeyType& key)
{
    HASH_INDEX_T h = probe(key);
    if(h == npos || table_[h] == nullptr) return nullptr;
    return &table_[h]->item;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
typename HashTable<K,V,Prober,Hash,KEqual>::HashItem*
HashTable<K,V,Prober,Hash,KEqual>::internalFind(const KeyType& key) const
{
    HASH_INDEX_T h = probe(key);
    if(h == npos || table_[h] == nullptr) return nullptr;
    return table_[h];
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::resize()
{
    if(mIndex_ + 1 >= sizeof(CAPACITIES)/sizeof(HASH_INDEX_T)) {
        throw std::logic_error("No more sizes");
    }

    std::vector<HashItem*> old = table_;

    mIndex_++;
    table_.clear();
    table_.resize(CAPACITIES[mIndex_], nullptr);

    currSize_ = 0;
    occupied_ = 0;

    
    double oldAlpha = resizeAlpha_;
    resizeAlpha_ = 1.0;

    for(size_t i = 0; i < old.size(); i++) {
        if(old[i] != nullptr) {
            if(!old[i]->deleted) {
                insert(old[i]->item);
            }
            delete old[i];
        }
    }

    resizeAlpha_ = oldAlpha;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
HASH_INDEX_T HashTable<K,V,Prober,Hash,KEqual>::probe(const KeyType& key) const
{
    HASH_INDEX_T h = hash_(key) % CAPACITIES[mIndex_];
    prober_.init(h, CAPACITIES[mIndex_], key);

    HASH_INDEX_T loc = prober_.next();
    totalProbes_++;

    while(loc != Prober::npos)
    {
        if(table_[loc] == nullptr) {
            return loc;
        }
        else if(!table_[loc]->deleted &&
                kequal_(table_[loc]->item.first, key))
        {
            return loc;
        }
        loc = prober_.next();
        totalProbes_++;
    }
    return npos;
}

template<typename K, typename V, typename Prober, typename Hash, typename KEqual>
void HashTable<K,V,Prober,Hash,KEqual>::reportAll(std::ostream& out) const
{
    for(HASH_INDEX_T i = 0; i < CAPACITIES[mIndex_]; ++i)
    {
        if(table_[i] != nullptr)
        {
            out << "Bucket " << i << ": "
                << table_[i]->item.first << " "
                << table_[i]->item.second << std::endl;
        }
    }
}

#endif
