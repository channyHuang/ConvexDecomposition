/* Copyright (c) 2011 Khaled Mamou (kmamou at gmail dot com)
All rights reserved.


 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

 3. The names of the contributors may not be used to endorse or promote products derived from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#pragma once
#ifndef VHACD_MANIFOLD_MESH_H
#define VHACD_MANIFOLD_MESH_H
#include "vhacdCircularList.h"
#include "vhacdSArray.h"
#include "vhacdVector.h"
namespace VHACD {
class TMMTriangle;
class TMMEdge;
class TMMesh;
class ICHull;

class DPoint
{
public:
                                                            DPoint(long name=0, Real dist=0, bool computed=false, bool distOnly=false)
                                                                :m_name(name),
                                                                 m_dist(dist),
                                                                 m_computed(computed),
                                                                 m_distOnly(distOnly){};
                                                            ~DPoint(){};
//private:
    long													m_name;
    Real                                                    m_dist;
    bool													m_computed;
    bool                                                    m_distOnly;
    /*friend class TMMTriangle;
    friend class TMMesh;
    friend class GraphVertex;
    friend class GraphEdge;
    friend class Graph;
    friend class ICHull;
    friend class HACD;*/
};

//!    Vertex data structure used in a triangular manifold mesh (TMM).
class TMMVertex {
public:
    void Initialize();
    TMMVertex(void);
    ~TMMVertex(void);

//private:
    Vec3<double> m_pos;
    int32_t m_name;
    size_t m_id;
    CircularListElement<TMMEdge>* m_duplicate; // pointer to incident cone edge (or NULL)
    bool m_onHull;
    bool m_tag;
    TMMVertex(const TMMVertex& rhs);
    friend class ICHull;
    friend class TMMesh;
    friend class TMMTriangle;
    friend class TMMEdge;
};

//!    Edge data structure used in a triangular manifold mesh (TMM).
class TMMEdge {
public:
    void Initialize();
    TMMEdge(void);
    ~TMMEdge(void);

//private:
    size_t m_id;
    CircularListElement<TMMTriangle>* m_triangles[2];
    CircularListElement<TMMVertex>* m_vertices[2];
    CircularListElement<TMMTriangle>* m_newFace;
    TMMEdge(const TMMEdge& rhs);
    friend class ICHull;
    friend class TMMTriangle;
    friend class TMMVertex;
    friend class TMMesh;
};

//!    Triangle data structure used in a triangular manifold mesh (TMM).
class TMMTriangle {
public:
    void Initialize();
    TMMTriangle(void);
    ~TMMTriangle(void);

//private:
    size_t m_id;
    CircularListElement<TMMEdge>* m_edges[3];
    CircularListElement<TMMVertex>* m_vertices[3];
    bool m_visible;

    SArray<long, SARRAY_DEFAULT_MIN_SIZE>				m_incidentPoints;

    TMMTriangle(const TMMTriangle& rhs);
    friend class ICHull;
    friend class TMMesh;
    friend class TMMVertex;
    friend class TMMEdge;
};
//!    triangular manifold mesh data structure.
class TMMesh {
public:
    //! Returns the number of vertices>
    inline size_t GetNVertices() const { return m_vertices.GetSize(); }
    //! Returns the number of edges
    inline size_t GetNEdges() const { return m_edges.GetSize(); }
    //! Returns the number of triangles
    inline size_t GetNTriangles() const { return m_triangles.GetSize(); }
    //! Returns the vertices circular list
    inline const CircularList<TMMVertex>& GetVertices() const { return m_vertices; }
    //! Returns the edges circular list
    inline const CircularList<TMMEdge>& GetEdges() const { return m_edges; }
    //! Returns the triangles circular list
    inline const CircularList<TMMTriangle>& GetTriangles() const { return m_triangles; }
    //! Returns the vertices circular list
    inline CircularList<TMMVertex>& GetVertices() { return m_vertices; }
    //! Returns the edges circular list
    inline CircularList<TMMEdge>& GetEdges() { return m_edges; }
    //! Returns the triangles circular list
    inline CircularList<TMMTriangle>& GetTriangles() { return m_triangles; }
    //! Add vertex to the mesh
    CircularListElement<TMMVertex>* AddVertex() { return m_vertices.Add(); }
    //! Add vertex to the mesh
    CircularListElement<TMMEdge>* AddEdge() { return m_edges.Add(); }
    //! Add vertex to the mesh
    CircularListElement<TMMTriangle>* AddTriangle() { return m_triangles.Add(); }
    //! Print mesh information
    void Print();
    //!
    void GetIFS(Vec3<double>* const points, Vec3<int32_t>* const triangles);
    //!
    void Clear();
    //!
    void Copy(TMMesh& mesh);
    //!
    bool CheckConsistancy();
    //!
    bool Normalize();
    //!
    bool Denormalize();
    //!    Constructor
    TMMesh();
    //! Destructor
    virtual ~TMMesh(void);

//private:
    CircularList<TMMVertex> m_vertices;
    CircularList<TMMEdge> m_edges;
    CircularList<TMMTriangle> m_triangles;

    // not defined
    TMMesh(const TMMesh& rhs);
    friend class ICHull;

};

long														IntersectRayTriangle( const Vec3<double> & P0, const Vec3<double> & dir,
                                                                                  const Vec3<double> & V0, const Vec3<double> & V1,
                                                                                  const Vec3<double> & V2, double &t);
}
#endif // VHACD_MANIFOLD_MESH_H
