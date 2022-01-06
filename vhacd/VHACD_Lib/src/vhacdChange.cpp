#include "vhacdChange.h"

namespace VHACD {


void VHACDChange::CreateGraph()
{
    // vertex to triangle adjacency information
    std::vector< std::set<long> >  vertexToTriangles;
    vertexToTriangles.resize(m_nPoints);
    for(size_t t = 0; t < m_nTriangles; ++t)
    {
        vertexToTriangles[m_triangles[t].X()].insert(static_cast<long>(t));
        vertexToTriangles[m_triangles[t].Y()].insert(static_cast<long>(t));
        vertexToTriangles[m_triangles[t].Z()].insert(static_cast<long>(t));
    }

    m_graph.Clear();
    m_graph.Allocate(m_nTriangles, 5 * m_nTriangles);
    unsigned long long tr1[3];
    unsigned long long tr2[3];
    long i1, j1, k1, i2, j2, k2;
    long t1, t2;
    for (size_t v = 0; v < m_nPoints; v++)
    {
        std::set<long>::const_iterator it1(vertexToTriangles[v].begin()), itEnd(vertexToTriangles[v].end());
        for(; it1 != itEnd; ++it1)
        {
            t1 = *it1;
            i1 = m_triangles[t1].X();
            j1 = m_triangles[t1].Y();
            k1 = m_triangles[t1].Z();
            tr1[0] = GetEdgeIndex(i1, j1);
            tr1[1] = GetEdgeIndex(j1, k1);
            tr1[2] = GetEdgeIndex(k1, i1);
            std::set<long>::const_iterator it2(it1);
            for(++it2; it2 != itEnd; ++it2)
            {
                t2 = *it2;
                i2 = m_triangles[t2].X();
                j2 = m_triangles[t2].Y();
                k2 = m_triangles[t2].Z();
                tr2[0] = GetEdgeIndex(i2, j2);
                tr2[1] = GetEdgeIndex(j2, k2);
                tr2[2] = GetEdgeIndex(k2, i2);
                int shared = 0;
                for(int i = 0; i < 3; ++i)
                {
                    for(int j = 0; j < 3; ++j)
                    {
                        if (tr1[i] == tr2[j])
                        {
                            shared++;
                        }
                    }
                }
                if (shared == 1) // two triangles are connected if they share exactly one edge
                {
                    m_graph.AddEdge(t1, t2);
                }
            }
        }
    }
    if (m_ccConnectDist >= 0.0)
    {
        m_graph.ExtractCCs();
        if (m_callBack)
        {
            char msg[1024];
            sprintf(msg, "nCC %lu\n", m_graph.m_nCCs);
            (*m_callBack)(msg, 0.0, 0.0,  m_graph.GetNVertices());

        }

        if (m_graph.m_nCCs > 1)
        {
            std::vector< std::set<long> > cc2V;
            cc2V.resize(m_graph.m_nCCs);
            long cc;
            for(size_t t = 0; t < m_nTriangles; ++t)
            {
                cc = m_graph.m_vertices[t].m_cc;
                cc2V[cc].insert(m_triangles[t].X());
                cc2V[cc].insert(m_triangles[t].Y());
                cc2V[cc].insert(m_triangles[t].Z());
            }

            for(size_t cc1 = 0; cc1 < m_graph.m_nCCs; ++cc1)
            {
                for(size_t cc2 = cc1+1; cc2 < m_graph.m_nCCs; ++cc2)
                {
                    std::set<long>::const_iterator itV1(cc2V[cc1].begin()), itVEnd1(cc2V[cc1].end());
                    for(; itV1 != itVEnd1; ++itV1)
                    {
                        double distC1C2 = std::numeric_limits<double>::max();
                        double dist;
                        t1 = -1;
                        t2 = -1;
                        std::set<long>::const_iterator itV2(cc2V[cc2].begin()), itVEnd2(cc2V[cc2].end());
                        for(; itV2 != itVEnd2; ++itV2)
                        {
                            dist = (m_points[*itV1] - m_points[*itV2]).GetNorm();
                            if (dist < distC1C2)
                            {
                                distC1C2 = dist;
                                t1 = *vertexToTriangles[*itV1].begin();

                                std::set<long>::const_iterator it2(vertexToTriangles[*itV2].begin()),
                                                               it2End(vertexToTriangles[*itV2].end());
                                t2 = -1;
                                for(; it2 != it2End; ++it2)
                                {
                                    if (*it2 != t1)
                                    {
                                        t2 = *it2;
                                        break;
                                    }
                                }
                            }
                        }
                        if (distC1C2 <= m_ccConnectDist && t1 >= 0 && t2 >= 0)
                        {
                            m_graph.AddEdge(t1, t2);
                        }
                    }
                }
            }
        }
    }
}

void VHACDChange::InitializeDualGraph()
{
    long i, j, k;
    Vec3<Real> u, v, w, normal;
    //delete [] m_normals;
    m_normals = new Vec3<Real>[m_nPoints];
    if (m_addFacesPoints)
    {
        //delete [] m_facePoints;
        //delete [] m_faceNormals;

        m_facePoints = new Vec3<Real>[m_nTriangles];
        m_faceNormals = new Vec3<Real>[m_nTriangles];
/*
        m_facePoints = new Vec3<Real>[4*m_nTriangles];
        m_faceNormals = new Vec3<Real>[4*m_nTriangles];
*/
    }
    memset(m_normals, 0, sizeof(Vec3<Real>) * m_nPoints);

    HACD::RaycastMesh rm;
    if (m_addExtraDistPoints)
    {
        rm.Initialize(m_nPoints, m_nTriangles, m_points, m_triangles);
        m_extraDistPoints = new Vec3<Real>[m_nTriangles];
        m_extraDistNormals = new Vec3<Real>[m_nTriangles];
    }
    double progressOld = -1.0;
    double progress = 0.0;
    char msg[1024];
    double ptgStep = 1.0;
    m_area = 0.0;
    for(unsigned long f = 0; f < m_nTriangles; f++)
    {
/*
        progress = f * 100.0 / m_nTriangles;
        if (fabs(progress-progressOld) > ptgStep && m_callBack)
        {
            sprintf(msg, "%3.2f %% \t \t \r", progress);
            (*m_callBack)(msg, progress, 0.0,  m_nTriangles);
            progressOld = progress;
        }
*/
        i = m_triangles[f].X();
        j = m_triangles[f].Y();
        k = m_triangles[f].Z();
        m_graph.m_vertices[f].m_distPoints.PushBack(DPoint(i, 0, false, false));
        m_graph.m_vertices[f].m_distPoints.PushBack(DPoint(j, 0, false, false));
        m_graph.m_vertices[f].m_distPoints.PushBack(DPoint(k, 0, false, false));

        HACD::ICHull  * ch = new HACD::ICHull(m_heapManager);
        m_graph.m_vertices[f].m_convexHull = ch;
        ch->AddPoint(m_points[i], i);
        ch->AddPoint(m_points[j], j);
        ch->AddPoint(m_points[k], k);
        ch->SetDistPoints(0);

        u = m_points[j] - m_points[i];
        v = m_points[k] - m_points[i];
        w = m_points[k] - m_points[j];
        normal = u ^ v;

        m_normals[i] += normal;
        m_normals[j] += normal;
        m_normals[k] += normal;

        m_graph.m_vertices[f].m_surf = normal.GetNorm();
        m_area += m_graph.m_vertices[f].m_surf;
        normal.Normalize();
        m_graph.m_vertices[f].m_boudaryEdges.Insert(GetEdgeIndex(i,j));
        m_graph.m_vertices[f].m_boudaryEdges.Insert(GetEdgeIndex(j,k));
        m_graph.m_vertices[f].m_boudaryEdges.Insert(GetEdgeIndex(k,i));
        if(m_addFacesPoints)
        {
            m_faceNormals[f] = normal;
            m_facePoints[f] = (m_points[i] + m_points[j] + m_points[k]) / 3.0;
            m_graph.m_vertices[f].m_distPoints.PushBack(DPoint(-static_cast<long>(f)-1, 0, false, true));
        }
    }

    // Now we have all the points in the KD tree, optimize the distance points insertion by running them in parallel
    // if possible.
    if (m_addExtraDistPoints)
    {/*
        if (m_callBack)
            (*m_callBack)("++ Also adding distance points\n", 0.0, 0.0, 0);
*/
        progressOld = -1.0;
        progress = 0.0;
        long completed = 0;

#ifdef THREAD_DIST_POINTS
#pragma omp parallel for
#endif

#ifdef CHANGE_SAMPLE
        Vec3<Real> vdir[3] = {m_points[i] - m_points[j], m_points[i] - m_points[k], m_points[j] - m_points[k]};
        for (int cdir = 0; cdir < 3; cdir++) {
            int steps = std::sqrt(vdir[cdir].X() * vdir[cdir].X() + vdir[cdir].Y() * vdir[cdir].Y() + vdir[cdir].Z() * vdir[cdir].Z());
            for (int sidx = 0; sidx < steps; sidx++) {
                Vec3<Real> seedPoint = (cdir == 2 ? m_points[k] : m_points[i]) + vdir[cdir] * sidx;
#endif
        for(long f = 0; f < (long)m_nTriangles; f++)
        {
#ifndef CHANGE_SAMPLE
            Vec3<Real> seedPoint((m_points[i] + m_points[j] + m_points[k]) / 3.0);
#endif
            Vec3<Real> hitPoint;
            Vec3<Real> hitNormal;
            normal = -normal;
            size_t faceIndex = m_nTriangles;

            HACD::Float dist;
            long hitTriangle;
            if (rm.Raycast(seedPoint, normal, hitTriangle,dist, hitPoint, hitNormal))
            {
                faceIndex = hitTriangle;
            }

            if (faceIndex < m_nTriangles )
            {
                m_extraDistPoints[f] = hitPoint;
                m_extraDistNormals[f] = hitNormal;
                m_graph.m_vertices[f].m_distPoints.PushBack(DPoint(m_nPoints+f, 0, false, true));
            }

            // Atomic update of the progress
            #ifdef THREAD_DIST_POINTS
            #pragma omp critical
            #endif
            {
                completed++;
/*
                progress = completed * 100.0 / m_nTriangles;
                if (fabs(progress-progressOld) > ptgStep && m_callBack)
                {
                    sprintf(msg, "%3.2f %% \t \t \r", progress);
                    (*m_callBack)(msg, progress, 0.0,  m_nTriangles);
                    progressOld = progress;
                }*/
            }
        }

#ifdef CHANGE_SAMPLE
            }
        }
#endif
    }

    for (size_t v = 0; v < m_nPoints; v++)
    {
        m_normals[v].Normalize();
    }
}

bool VHACDChange::InitializePriorityQueue()
{
    m_pqueue.reserve(m_graph.m_nE + 100);
    for (size_t e=0; e < m_graph.m_nE; ++e)
    {
        ComputeEdgeCost(static_cast<long>(e));
        m_pqueue.push(GraphEdgePriorityQueue(static_cast<long>(e), m_graph.m_edges[e].m_error));
    }
    return true;
}

void VHACDChange::ComputeEdgeCost(size_t e)
{
    GraphEdge & gE = m_graph.m_edges[e];
    long v1 = gE.m_v1;
    long v2 = gE.m_v2;

    if (m_graph.m_vertices[v2].m_ancestors.size()>m_graph.m_vertices[v1].m_ancestors.size())
    {
        gE.m_v1 = v2;
        gE.m_v2 = v1;
        std::swap(v1, v2);
    }
    GraphVertex & gV1 = m_graph.m_vertices[v1];
    GraphVertex & gV2 = m_graph.m_vertices[v2];
#ifdef HACD_DEBUG
    if (v1 == 308 && v2==276)
    {
        gV1.m_convexHull->m_mesh.Save("debug1.wrl");
        gV2.m_convexHull->m_mesh.Save("debug2.wrl");
    }

#endif

    // create the edge's convex-hull
    HACD::ICHull  * ch = new HACD::ICHull(m_heapManager);
    (*ch) = (*gV1.m_convexHull);
    // update distPoints
#ifdef HACD_PRECOMPUTE_CHULLS
    delete gE.m_convexHull;
    gE.m_convexHull = 0;
#endif
    std::map<long, DPoint> distPoints;
    for(size_t p = 0; p < gV1.m_distPoints.Size(); ++p)
    {
        distPoints[gV1.m_distPoints[p].m_name] = gV1.m_distPoints[p];
    }

    std::map<long, DPoint>::iterator itDP1;
    for(size_t p = 0; p < gV2.m_distPoints.Size(); ++p)
    {
        const DPoint & point =  gV2.m_distPoints[p];
        itDP1 = distPoints.find(point.m_name);
        if (itDP1 == distPoints.end())
        {
            DPoint newPoint(point.m_name, 0, false, point.m_distOnly);
            distPoints.insert(std::pair<long, DPoint>(point.m_name, newPoint));
            if ( !point.m_distOnly )
            {
                ch->AddPoint(m_points[point.m_name], point.m_name);
            }
        }
        else
        {
            if ( (itDP1->second).m_distOnly && !point.m_distOnly)
            {
                (itDP1->second).m_distOnly = false;
                ch->AddPoint(m_points[point.m_name], point.m_name);
            }
        }
    }

    ch->SetDistPoints(&distPoints);
    // create the convex-hull
    while (ch->Process() == HACD::ICHullErrorInconsistent)		// if we face problems when constructing the visual-hull. really ugly!!!!
    {
//			if (m_callBack) (*m_callBack)("\t Problem with convex-hull construction [HACD::ComputeEdgeCost]\n", 0.0, 0.0, 0);
        HACD::ICHull  * chOld = ch;
        ch = new HACD::ICHull(m_heapManager);
        CircularList<TMMVertex> & verticesCH = chOld->GetMesh().m_vertices;
        size_t nV = verticesCH.GetSize();
        long ptIndex = 0;
        verticesCH.Next();
        // add noise to avoid the problem
        ptIndex = verticesCH.GetHead()->GetData().m_name;
        ch->AddPoint(m_points[ptIndex]+ m_scale * 0.0001 * Vec3<Real>(rand() % 10 - 5, rand() % 10 - 5, rand() % 10 - 5), ptIndex);
        for(size_t v = 1; v < nV; ++v)
        {
            ptIndex = verticesCH.GetHead()->GetData().m_name;
            ch->AddPoint(m_points[ptIndex], ptIndex);
            verticesCH.Next();
        }
        delete chOld;
    }
#ifdef HACD_DEBUG
    if (v1 == 438 && v2==468)
    {
        const long nPoints = static_cast<long>(m_nPoints);
        std::map<long, DPoint>::iterator itDP(distPoints.begin());
        std::map<long, DPoint>::iterator itDPEnd(distPoints.end());
        for(; itDP != itDPEnd; ++itDP)
        {

            if (itDP->first >= nPoints)
            {
                long pt = itDP->first - nPoints;
                ch->AddPoint(m_extraDistPoints[pt], itDP->first);
            }
            else if (itDP->first >= 0)
            {
                long pt = itDP->first;
                ch->AddPoint(m_points[pt], itDP->first);
            }
            else
            {
                long pt = -itDP->first-1;
                ch->AddPoint(m_facePoints[pt], itDP->first);
                ch->AddPoint(m_facePoints[pt] + 10.0 * m_faceNormals[pt] , itDP->first);
            }
        }
        printf("-***->\n");

        ch->m_mesh.Save("debug.wrl");
    }
#endif
    double surf = gV1.m_surf + gV2.m_surf;
    double concavity = 0.0;
    double surfCH = ch->ComputeArea() / 2.0;
    double volumeCH = ch->ComputeVolume();
    double vol2Surf = volumeCH / surfCH;
    double concavity_flat = sqrt(fabs(surfCH-surf));
    double weightFlat = std::max(0.0, 1.0 - pow(- vol2Surf * 100.0 / (m_scale * m_flatRegionThreshold), 2.0));

//		concavity_flat *= std::max(exp(- vol2Surf * 100.0 / (m_scale * m_flatRegionThreshold)) - exp(-1.0), 0.0);
    concavity_flat *= weightFlat;
    if(!ch->IsFlat())
    {
        concavity = Concavity(*ch, distPoints);
#ifdef CHANGLE_CONCAVITY
        concavity = (gV1.m_convexHull->ComputeArea() + gV2.m_convexHull->ComputeArea() - ch->ComputeArea());
        //concavity = std::sqrt(concavity);
        concavity *= concavity;
#endif
    }
    concavity += concavity_flat;
#ifdef HACD_PRECOMPUTE_CHULLS
    gE.m_convexHull = ch;
#else
    delete ch;
#endif

    // compute boudary edges
    double perimeter = 0.0;
    if (m_alpha > 0.0)
    {
        std::set<unsigned long long> boudaryEdges1;
        for(size_t edV1 = 0; edV1 < gV1.m_boudaryEdges.Size(); ++edV1)
        {
            boudaryEdges1.insert(gV1.m_boudaryEdges[edV1]);
        }
        std::set<unsigned long long> boudaryEdges2;
        for(size_t edV2 = 0; edV2 < gV2.m_boudaryEdges.Size(); ++edV2)
        {
            boudaryEdges2.insert(gV2.m_boudaryEdges[edV2]);
        }
        std::set<unsigned long long> boudaryEdges;
        std::set_symmetric_difference (boudaryEdges1.begin(),
                                       boudaryEdges1.end(),
                                       boudaryEdges2.begin(),
                                       boudaryEdges2.end(),
                                       std::inserter( boudaryEdges, boudaryEdges.begin() ) );

        std::set<unsigned long long>::const_iterator itBE(boudaryEdges.begin());
        std::set<unsigned long long>::const_iterator itBEEnd(boudaryEdges.end());
        for(; itBE != itBEEnd; ++itBE)
        {
                perimeter += (m_points[static_cast<long>((*itBE) >> 32)] -
                               m_points[static_cast<long>((*itBE) & 0xFFFFFFFFULL)]).GetNorm();
        }
    }
    double ratio   = perimeter * perimeter / (4.0 * FLOAT_MATH::FM_PI * surf);
    gE.m_concavity = concavity;                     // cluster's concavity
    double volume  = volumeCH/pow(m_scale, 3.0);	// cluster's volume
    gE.m_error     = static_cast<Real>(concavity +  m_alpha * (1.0 - weightFlat) * ratio + m_beta * volume + m_gamma * static_cast<double>(distPoints.size()) / m_nPoints);	// cluster's priority
}

double VHACDChange::Concavity(HACD::ICHull & ch, std::map<long, DPoint> & distPoints)
{
    double concavity = 0.0;
    double distance = 0.0;
    std::map<long, DPoint>::iterator itDP(distPoints.begin());
    std::map<long, DPoint>::iterator itDPEnd(distPoints.end());
    long pt;
    const long nPoints = static_cast<long>(m_nPoints);
    const double eps = -0.001;
    for(; itDP != itDPEnd; ++itDP)
    {
        if (!(itDP->second).m_computed)
        {
            if (itDP->first >= nPoints)
            {
                pt = itDP->first - nPoints;
                if ( ch.IsInside(m_extraDistPoints[pt], eps))
                {
                    distance = ch.ComputeDistance(itDP->first, m_extraDistPoints[pt], m_extraDistNormals[pt], (itDP->second).m_computed, true);
                }
                else
                {
                    distance = 0.0;
                }
            }
            else if (itDP->first >= 0)
            {
                pt = itDP->first;
                distance = ch.ComputeDistance(itDP->first, m_points[pt], m_normals[pt], (itDP->second).m_computed, true);
            }
            else
            {
                pt = -itDP->first-1;
                distance = ch.ComputeDistance(itDP->first, m_facePoints[pt], m_faceNormals[pt], (itDP->second).m_computed, true);
            }
        }
        else
        {
            distance = (itDP->second).m_dist;
        }

        if (concavity < distance)
        {
            concavity = distance;
        }
    }
    return concavity;
}

void VHACDChange::GetClipPlanes(SArray<Plane> &planes) {
    auto tmpQuque = m_pqueue;
    while (!tmpQuque.empty()) {
        GraphEdgePriorityQueue currentEdge(0,0.0);
        currentEdge = tmpQuque.top();
        tmpQuque.pop();

        long v1 = m_graph.m_edges[currentEdge.m_name].m_v1;
        long v2 = m_graph.m_edges[currentEdge.m_name].m_v2;

        {
            long index[3] = { m_triangles[v1].X(), m_triangles[v1].Y(), m_triangles[v1].Z() };
            auto dir1 = m_points[index[1]] - m_points[index[0]];
            auto dir2 = m_points[index[2]] - m_points[index[0]];
            auto normal1 = (dir1 ^ dir2);

            Plane plane;
            double param[3];
            double p1[3] = { m_points[index[0]].X(), m_points[index[0]].Y(), m_points[index[0]].Z() };
            double p2[3] = { m_points[index[1]].X(), m_points[index[1]].Y(), m_points[index[1]].Z() };
            double p3[3] = { m_points[index[2]].X(), m_points[index[2]].Y(), m_points[index[2]].Z() };
            plane.m_d = FLOAT_MATH::fm_computePlane(p1, p2, p3, param);
            plane.m_a = param[0];
            plane.m_b = param[1];
            plane.m_c = param[2];

            planes.PushBack(plane);
        }

        {
            long index[3] = { m_triangles[v2].X(), m_triangles[v2].Y(), m_triangles[v2].Z() };
            auto dir1 = m_points[index[1]] - m_points[index[0]];
            auto dir2 = m_points[index[2]] - m_points[index[0]];
            auto normal1 = (dir1 ^ dir2);

            Plane plane;
            double param[3];
            double p1[3] = { m_points[index[0]].X(), m_points[index[0]].Y(), m_points[index[0]].Z() };
            double p2[3] = { m_points[index[1]].X(), m_points[index[1]].Y(), m_points[index[1]].Z() };
            double p3[3] = { m_points[index[2]].X(), m_points[index[2]].Y(), m_points[index[2]].Z() };
            plane.m_d = FLOAT_MATH::fm_computePlane(p1, p2, p3, param);
            plane.m_a = param[0];
            plane.m_b = param[1];
            plane.m_c = param[2];

            planes.PushBack(plane);
        }
    }
}

}
