#ifndef VHACDCHANGE_H
#define VHACDCHANGE_H

#include "hacdGraph.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#if _OPENMP
#include <omp.h>
#endif // _OPENMP

#include <vector>
#include <set>
#include <queue>

//#include "vhacdRaycastMesh.h"
#include "hacdICHull.h"
#include "FloatMath.h"
#include "hacdRaycastMesh.h"
#include "vhacdMesh.h"
#include "hacdMicroAllocator.h"
namespace VHACD {

template<class _Ty, class _Container = std::vector<_Ty>, class _Pr = std::less<typename _Container::value_type> >
class reservable_priority_queue: public std::priority_queue<_Ty, _Container, _Pr>
{
    typedef typename std::priority_queue<_Ty, _Container, _Pr>::size_type size_type;
public:
                                                reservable_priority_queue(size_type capacity = 0) { reserve(capacity); };
    void										reserve(size_type capacity) { this->c.reserve(capacity); }
    size_type									capacity() const { return this->c.capacity(); }
};

class GraphEdgePriorityQueue
{
    public:
        //! Constructor
        //! @param name edge's id
        //! @param priority edge's priority
                                                GraphEdgePriorityQueue(long name, Real priority)
                                                {
                                                    m_name = name;
                                                    m_priority = priority;
                                                }
        //! Destructor
                                                ~GraphEdgePriorityQueue(void){}
//    private:
        long									m_name;						//!< edge name
        Real                                    m_priority;					//!< priority
    //! Operator < for GraphEdgePQ
    friend bool                                 operator<(const GraphEdgePriorityQueue & lhs, const GraphEdgePriorityQueue & rhs);
    //! Operator > for GraphEdgePQ
    friend bool                                 operator>(const GraphEdgePriorityQueue & lhs, const GraphEdgePriorityQueue & rhs);
    //friend class VHACD;
};

inline bool										operator<(const GraphEdgePriorityQueue & lhs, const GraphEdgePriorityQueue & rhs)
                                                {
                                                    return (lhs.m_priority<rhs.m_priority);
                                                }
inline bool										operator>(const GraphEdgePriorityQueue & lhs, const GraphEdgePriorityQueue & rhs)
                                                {
                                                    return lhs.m_priority>rhs.m_priority;
                                                }

class VHACDChange
{
public:
    typedef void (*CallBackFunction)(const char *, double, double, size_t);


    VHACDChange() {};
    ~VHACDChange() {};

    void CreateGraph();
    void InitializeDualGraph();
    bool InitializePriorityQueue();
    void ComputeEdgeCost(size_t e);
    double Concavity(HACD::ICHull & ch, std::map<long, DPoint> & distPoints);
    void GetClipPlanes(SArray<Plane> &planes);

    static unsigned long long					GetEdgeIndex(unsigned long long a, unsigned long long b)
                                                {
                                                    if (a > b) return (a << 32) + b;
                                                    else	   return (b << 32) + a;
                                                }

    //! Gives the targeted number of triangles of the decimated mesh
    //! @return targeted number of triangles of the decimated mesh
    size_t										GetTargetNTrianglesDecimatedMesh() const { return m_targetNTrianglesDecimatedMesh;}
    //! Sets the targeted number of triangles of the decimated mesh
    //! @param targeted number of triangles of the decimated mesh
    void										SetNTargetTrianglesDecimatedMesh(size_t  targetNTrianglesDecimatedMesh) { m_targetNTrianglesDecimatedMesh = targetNTrianglesDecimatedMesh;}
    //! Gives the triangles partitionas an array of size m_nTriangles where the i-th element specifies the cluster to which belong the i-th triangle
    //! @return triangles partition
    const long * const							GetPartition() const { return m_partition;}
    //! Sets the scale factor
    //! @param scale scale factor
    void										SetScaleFactor(double  scale) { m_scale = scale;}
    //! Gives the scale factor
    //! @return scale factor
    const double								GetScaleFactor() const { return m_scale;}
    //! Sets the threshold to detect small clusters. The threshold is expressed as a percentage of the total mesh surface (default 0.25%)
    //! @param smallClusterThreshold threshold to detect small clusters
    void										SetSmallClusterThreshold(double  smallClusterThreshold) { m_smallClusterThreshold = smallClusterThreshold;}
    //! Gives the threshold to detect small clusters. The threshold is expressed as a percentage of the total mesh surface (default 0.25%)
    //! @return threshold to detect small clusters
    const double								GetSmallClusterThreshold() const { return m_smallClusterThreshold;}
    //! Sets the call-back function
    //! @param callBack pointer to the call-back function
    void										SetCallBack(CallBackFunction  callBack) { m_callBack = callBack;}
    //! Gives the call-back function
    //! @return pointer to the call-back function
    const CallBackFunction                      GetCallBack() const { return m_callBack;}

    //! Specifies whether faces points should be added when computing the concavity
    //! @param addFacesPoints true = faces points should be added
    void										SetAddFacesPoints(bool  addFacesPoints) { m_addFacesPoints = addFacesPoints;}
    //! Specifies wheter faces points should be added when computing the concavity
    //! @return true = faces points should be added
    const bool									GetAddFacesPoints() const { return m_addFacesPoints;}
    //! Specifies whether extra points should be added when computing the concavity
    //! @param addExteraDistPoints true = extra points should be added
    void										SetAddExtraDistPoints(bool  addExtraDistPoints) { m_addExtraDistPoints = addExtraDistPoints;}
    //! Specifies wheter extra points should be added when computing the concavity
    //! @return true = extra points should be added
    const bool									GetAddExtraDistPoints() const { return m_addExtraDistPoints;}
    //! Sets the points of the input mesh (Remark: the input points will be scaled and shifted. Use DenormalizeData() to invert those operations)
    //! @param points pointer to the input points
    void										SetPoints(Vec3<Real>  * points) { m_points = points;}
    //! Gives the points of the input mesh (Remark: the input points will be scaled and shifted. Use DenormalizeData() to invert those operations)
    //! @return pointer to the input points
    const Vec3<Real> *                          GetPoints() const { return m_points;}
    //! Gives the points of the decimated mesh
    //! @return pointer to the decimated mesh points
    const Vec3<Real> *                          GetDecimatedPoints() const { return m_pointsDecimated;}
    //! Gives the triangles in the decimated mesh
    //! @return pointer to the decimated mesh triangles
    const Vec3<long>   *			            GetDecimatedTriangles() const { return m_trianglesDecimated;}
    //! Gives the number of points in the decimated mesh.
    //! @return number of points the decimated mesh mesh
    const size_t								GetNDecimatedPoints() const { return m_nPointsDecimated;}
    //! Gives the number of triangles in the decimated mesh.
    //! @return number of triangles the decimated mesh
    const size_t								GetNDecimatedTriangles() const { return m_nTrianglesDecimated;}
    //! Sets the triangles of the input mesh.
    //! @param triangles points pointer to the input points
    void										SetTriangles(Vec3<long>  * triangles) { m_triangles = triangles;}
    //! Gives the triangles in the input mesh
    //! @return pointer to the input triangles
    const Vec3<long>   *			            GetTriangles() const { return m_triangles;}
    //! Sets the number of points in the input mesh.
    //! @param nPoints number of points the input mesh
    void										SetNPoints(size_t nPoints) { m_nPoints = nPoints;}
    //! Gives the number of points in the input mesh.
    //! @return number of points the input mesh
    const size_t								GetNPoints() const { return m_nPoints;}
    //! Sets the number of triangles in the input mesh.
    //! @param nTriangles number of triangles in the input mesh
    void										SetNTriangles(size_t nTriangles) { m_nTriangles = nTriangles;}
    //! Gives the number of triangles in the input mesh.
    //! @return number of triangles the input mesh
    const size_t								GetNTriangles() const { return m_nTriangles;}
    //! Sets the minimum number of clusters to be generated.
    //! @param nClusters minimum number of clusters
    void										SetNClusters(size_t nClusters) { m_nMinClusters = nClusters;}
    //! Gives the number of generated clusters.
    //! @return number of generated clusters
    const size_t								GetNClusters() const { return m_nClusters;}
    //! Sets the maximum allowed concavity.
    //! @param concavity maximum concavity
    void										SetConcavity(double concavity) { m_concavity = concavity;}
    //! Gives the maximum allowed concavity.
    //! @return maximum concavity
    double                                      GetConcavity() const { return m_concavity;}
    //! Sets the maximum allowed distance to get CCs connected.
    //! @param concavity maximum distance to get CCs connected
    void										SetConnectDist(double ccConnectDist) { m_ccConnectDist = ccConnectDist;}
    //! Gives the maximum allowed distance to get CCs connected.
    //! @return maximum distance to get CCs connected
    double                                      GetConnectDist() const { return m_ccConnectDist;}
    //! Sets the volume weight.
    //! @param beta volume weight
    void										SetVolumeWeight(double beta) { m_beta = beta;}
    //! Gives the volume weight.
    //! @return volume weight
    double                                      GetVolumeWeight() const { return m_beta;}
    //! Sets the compacity weight (i.e. parameter alpha in ftp://ftp.elet.polimi.it/users/Stefano.Tubaro/ICIP_USB_Proceedings_v2/pdfs/0003501.pdf).
    //! @param alpha compacity weight
    void										SetCompacityWeight(double alpha) { m_alpha = alpha;}
    //! Gives the compacity weight (i.e. parameter alpha in ftp://ftp.elet.polimi.it/users/Stefano.Tubaro/ICIP_USB_Proceedings_v2/pdfs/0003501.pdf).
    //! @return compacity weight
    double                                      GetCompacityWeight() const { return m_alpha;}
    //! Sets the maximum number of vertices for each generated convex-hull.
    //! @param nVerticesPerCH maximum # vertices per CH
    void										SetNVerticesPerCH(size_t nVerticesPerCH) { m_nVerticesPerCH = nVerticesPerCH;}
    //! Gives the maximum number of vertices for each generated convex-hull.
    //! @return maximum # vertices per CH
    const size_t								GetNVerticesPerCH() const { return m_nVerticesPerCH;}
    //! Gives the number of vertices for the cluster number numCH.
    //! @return number of vertices

private:
    Vec3<long> *								m_trianglesDecimated;		//>! pointer the triangles array
    Vec3<Real> *                                m_pointsDecimated;			//>! pointer the points array
    Vec3<long> *								m_triangles;				//>! pointer the triangles array
    Vec3<Real> *                                m_points;					//>! pointer the points array
    Vec3<Real> *                                m_facePoints;               //>! pointer to the faces points array
    Vec3<Real> *                                m_faceNormals;              //>! pointer to the faces normals array
    Vec3<Real> *                                m_normals;					//>! pointer the normals array
    Vec3<Real> *                                m_extraDistPoints;          //>! pointer to the faces points array
    Vec3<Real> *                                m_extraDistNormals;         //>! pointer to the faces normals array
    size_t										m_nTrianglesDecimated;		//>! number of triangles in the original mesh
    size_t										m_nPointsDecimated;			//>! number of vertices in the original mesh
    size_t										m_nTriangles;				//>! number of triangles in the original mesh
    size_t										m_nPoints;					//>! number of vertices in the original mesh
    size_t										m_nClusters;				//>! number of clusters
    size_t										m_nMinClusters;				//>! minimum number of clusters
    double										m_ccConnectDist;			//>! maximum allowed distance to connect CCs
    double										m_concavity;				//>! maximum concavity
    double										m_alpha;					//>! compacity weigth
    double                                      m_beta;                     //>! volume weigth
    double										m_gamma;					//>! computation cost
    double										m_diag;						//>! length of the BB diagonal
    double										m_scale;					//>! scale factor used for NormalizeData() and DenormalizeData()
    double										m_flatRegionThreshold;		//>! threshhold to control the contirbution of flat regions concavity (default 1% of m_scale)
    double										m_smallClusterThreshold;	//>! threshhold to detect small clusters (default 0.25% of the total mesh surface)
    double										m_area;						//>! surface area
    Vec3<Real>                                  m_barycenter;				//>! barycenter of the mesh
    std::vector< long >                         m_cVertices;				//>! array of vertices each belonging to a different cluster
    HACD::ICHull *                                    m_convexHulls;				//>! convex-hulls associated with the final HACD clusters
    HACD::HeapManager *                               m_heapManager;              //>! Heap Manager
    Graph										m_graph;					//>! simplification graph
    size_t                                      m_nVerticesPerCH;			//>! maximum number of vertices per convex-hull

    reservable_priority_queue<GraphEdgePriorityQueue,
        std::vector<GraphEdgePriorityQueue>,
        std::greater<std::vector<GraphEdgePriorityQueue>::value_type> > m_pqueue;		//!> priority queue
//                                                VHACD(const VHACD & rhs);

    CallBackFunction							m_callBack;					//>! call-back function
    long *										m_partition;				//>! array of size m_nTriangles where the i-th element specifies the cluster to which belong the i-th triangle
    size_t										m_targetNTrianglesDecimatedMesh; //>! specifies the target number of triangles in the decimated mesh. If set to 0 no decimation is applied.
    //HeapManager *                               m_heapManager;              //>! Heap Manager
    bool                                        m_addFacesPoints;           //>! specifies whether to add faces points or not
    bool                                        m_addExtraDistPoints;       //>! specifies whether to add extra points for concave shapes or not

};

}
#endif
