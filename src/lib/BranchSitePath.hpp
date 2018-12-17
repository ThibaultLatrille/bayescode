#pragma once

#include <map>
#include <string>

class SuffStat;
class PoissonSuffStat;
class PathSuffStat;
class PolySuffStat;
class PathSuffStatBidimArray;
class SubMatrix;

/**
 * \brief The building block for substitution histories, or paths
 * (BranchSitePath class)
 *
 * A substitution history over some period of time, and with n substitution
 * events in total, is encoded as doubly-linked chains of n+1 Plink objects (see
 * BranchSitePath for an example). Each Plink encodes the current state and the
 * relative waiting time (relative to total time t) until either the next event
 * or the endpoint. Each Plink points to the previous Plink (Plink* prev) or to
 * the next (Plink* next). The first Plink has prev=null and the last Pling has
 * next=null;
 *
 * Thus, for instance
 * the following history, over a total time of t=10 units: A---C--T-----
 * (starting from A, waiting 3 time units, then making a substitution toward C,
 * then waiting 2 time units, making a substitution toward T, then waiting 10
 * time units and then stopping) would be encoded by a chain of three Plink
 * objects: (A,0.3) -> (C,0.2) -> (T,0.5).
 */

class Plink {
    friend class BranchSitePath;

  public:
    //! default constructor (empty)
    Plink();

    //! constructor specifying the current state and the relative time until next
    //! event
    Plink(int instate, double inrel_time);
    ~Plink();

    Plink *Prev();
    Plink *Next();

    const Plink *Prev() const;
    const Plink *Next() const;

    bool IsFirst() const;
    bool IsLast() const;

    void Splice();
    void Insert(Plink *link);

    void SetState(int instate);
    int GetState() const;

    void SetRelativeTime(double inrel_time);
    double GetRelativeTime() const;

  private:
    Plink *next;
    Plink *prev;

    int state;
    double rel_time;
};

/**
 * \brief A doubly-linked structure specifying the detailed substitution history
 * over a branch, for a given site
 *
 * A substitution history over some period of time, and with n substitution
 * events in total, is encoded as doubly-linked chains of n+1 Plink objects.
 * Thus, for instance
 * the following history, over a total time of t=10 units: A---C--T-----
 * (starting from A, waiting 3 time units, then making a substitution toward C,
 * then waiting 2 time units, making a substitution toward T, then waiting 10
 * time units and then stopping) would be encoded by a chain of three Plink
 * objects: (A,0.3) -> (C,0.2) -> (T,0.5).
 */

class BranchSitePath {
  public:
    //
    //! default constructor
    BranchSitePath();

    //! constructor for a minimal history: starting in given state, and then no
    //! event along entire branch
    BranchSitePath(int state);
    virtual ~BranchSitePath();

    //! const access to first Plink (at time 0)
    const Plink *Init() const;

    //! const access to last Plink (at the time of the last substitution event --
    //! can be the same as the first Plink if no event occured)
    const Plink *Last() const;

    //! non-const access to first Plink
    Plink *Init();
    //! non-const access to last Plink
    Plink *Last();

    //! return total number of substitution events
    int GetNsub() const;

    //! return initial state
    int GetInitState() const;

    //! return final state
    int GetFinalState() const;

    //! give the relative time for the event encoded by given Plink (0 if
    //! link==init)
    double GetRelativeTime(const Plink *link) const { return link->GetRelativeTime(); }

    //! \brief push up the sufficient statistics for this substitution history, as
    //! a function of the effective branch length, into the PoissonSuffStat given
    //! as first argument
    //!
    //! The sufficient statistic is: the total number of substitution events and
    //! the total effective rate over the whole branch. Note that this suff stat
    //! depends on the substitution process, such as specified by the SubMatrix.
    void AddLengthSuffStat(PoissonSuffStat &suffstat, double factor, const SubMatrix &mat) const;

    //! \brief push up the sufficient statistics for this substitution history, as
    //! a function of the substitution matrix, into the PathSuffStat given as
    //! first argument
    //!
    //! The sufficient statistic is: the total number of events from i->j for any
    //! (i,j) pair, and the total waiting time in state i, for each i. The factor
    //! given as second argument acts as a scaling factor for the waiting times.
    void AddPathSuffStat(PathSuffStat &suffstat, double factor) const;

    //! delete the current substitution history and create a new one starting with
    //! given initial state
    void Reset(int state);

    //! append a new event, after relative time reltimelength, and leading to new
    //! state instate
    void Append(int instate, double reltimelength);

    /*
    // reset the backup
    void BKReset(int state);
    // append a new event to the backup
    void BKAppend(int instate, double reltimelength);
    // backup current path
    void BackupPath();
    // restore current path, based on backup
    void RestorePath();
    */

  protected:
    Plink *init;
    Plink *last;
    int nsub;

    /*
    Plink *bkinit;
    Plink *bklast;
    int bknsub;
    */
};

//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//	* Inline definitions
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------

//-------------------------------------------------------------------------
//	* Plink
//-------------------------------------------------------------------------

inline Plink::Plink() : next(nullptr), prev(nullptr), state(0), rel_time(0) {}
inline Plink::Plink(int instate, double inrel_time)
    : next(nullptr), prev(nullptr), state(instate), rel_time(inrel_time) {}
inline Plink::~Plink() { Splice(); }

inline Plink *Plink::Prev() { return prev; }
inline Plink *Plink::Next() { return next; }

inline const Plink *Plink::Prev() const { return prev; }
inline const Plink *Plink::Next() const { return next; }

inline bool Plink::IsFirst() const { return prev == nullptr; }
inline bool Plink::IsLast() const { return next == nullptr; }

inline void Plink::Insert(Plink *link) {
    link->next = next;
    if (next != nullptr) { next->prev = link; }
    link->prev = this;
    next = link;
}

inline void Plink::Splice() {
    if (prev != nullptr) { prev->next = next; }
    if (next != nullptr) { next->prev = prev; }
    prev = next = nullptr;
}

inline void Plink::SetState(int instate) { state = instate; }
inline void Plink::SetRelativeTime(double inrel_time) { rel_time = inrel_time; }
inline double Plink::GetRelativeTime() const { return rel_time; }
inline int Plink::GetState() const { return state; }

inline Plink *BranchSitePath::Init() { return init; }
inline Plink *BranchSitePath::Last() { return last; }

inline const Plink *BranchSitePath::Init() const { return init; }
inline const Plink *BranchSitePath::Last() const { return last; }

//-------------------------------------------------------------------------
//	* BranchSitePath
//-------------------------------------------------------------------------

inline void BranchSitePath::Append(int instate, double reltimelength) {
    last->SetRelativeTime(reltimelength);
    auto link = new Plink(instate, 0);
    last->Insert(link);
    last = link;
    nsub++;
}

/*
inline void BranchSitePath::BKAppend(int instate, double reltimelength) {
    bklast->SetRelativeTime(reltimelength);
    auto link = new Plink(instate, 0);
    bklast->Insert(link);
    bklast = link;
    bknsub++;
}
*/

inline int BranchSitePath::GetNsub() const { return nsub; }

inline int BranchSitePath::GetInitState() const { return init->GetState(); }
inline int BranchSitePath::GetFinalState() const { return last->GetState(); }
