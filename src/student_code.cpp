/*
 * Student solution for CMU 15-462 Project 2 (MeshEdit)
 *
 * Implemented by ____ on ____.
 *
 */

#include "student_code.h"
#include "mutablePriorityQueue.h"

namespace CMU462
{
   VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
   {
      // TODO This method should split the given edge and return an iterator to the newly inserted vertex.
      // TODO The halfedge of this vertex should point along the edge that was split, rather than the new edges.

      
      // add three edges
      // recaculate the 4 part of the m esh

      // ## calculate the pos of new vertex
      // get half edge
      HalfedgeIter h = e0->halfedge();
      if (h->isBoundary()) {
         return VertexIter();
      }

      // get related two vertex
      VertexIter x0 = h->vertex();
      VertexIter x1 = h->twin()->vertex();

      HalfedgeIter first_h = h->next();
      HalfedgeIter first_h_twin = h->twin()->next();
      

      // change vertex and vertex degree
         // vertexDegree[h->vertex()]--;
         // vertexDegree[h_twin->vertex()]--;

      // h->vertex()->halfedge() = first_h_twin;
      // h->twin()->vertex()->halfedge() = first_h;

      // calculate new vertex position
      VertexIter midpoint = newVertex();
      midpoint->position = 0.5 * (x0->position + x1->position);

      // ## add new edges and related halfedges

      std::vector<HalfedgeIter> adjacent_halfedges;
      adjacent_halfedges.empty();
      adjacent_halfedges.push_back(h->next());
      adjacent_halfedges.push_back(h->next()->next());
      adjacent_halfedges.push_back(h->twin()->next());
      adjacent_halfedges.push_back(h->twin()->next()->next());

      std::vector<VertexIter> adjacent_vertices;
      adjacent_vertices.empty();
      for (int i = 0; i < adjacent_halfedges.size(); i++) 
         adjacent_vertices.push_back(adjacent_halfedges[i]->vertex());
      
      std::vector<FaceIter> related_faces;      
      related_faces.push_back(h->face());
      related_faces.push_back(h->twin()->face());
      for (int i = 0; i < 2; i++) related_faces.push_back(newFace());

      std::vector<EdgeIter> related_edges;
      related_edges.empty();
      related_edges.push_back(e0);
      for (int i = 0; i < 3; i++) related_edges.push_back(newEdge());

      std::vector<HalfedgeIter> related_halfedges;
      related_halfedges.empty();
      related_halfedges.push_back(h);
      related_halfedges.push_back(h->twin());
      for (int i = 0; i < 6; i++) related_halfedges.push_back(newHalfedge());

      // ## new face
      // FaceIter midpoint_face = newFace();
      // midpoint_face->halfedge() = related_halfedges[1];
      midpoint->halfedge() = related_halfedges[0];

      for (int i = 0; i < 4; i++) {
         related_halfedges[i * 2]->vertex() = midpoint;
         related_halfedges[i * 2 + 1]->vertex() = adjacent_vertices[i];

         related_halfedges[i * 2]->twin() = related_halfedges[i * 2 + 1];
         related_halfedges[i * 2 + 1]->twin() = related_halfedges[i * 2];

         related_halfedges[i * 2]->edge() = related_edges[i];
         related_halfedges[i * 2 + 1]->edge() = related_edges[i];

         related_halfedges[i * 2]->next() = adjacent_halfedges[i];
         related_halfedges[i * 2 + 1]->next() = related_halfedges[((i + 3) % 4) * 2];
         adjacent_halfedges[(i + 3) % 4]->next() = related_halfedges[i * 2 + 1];

         related_faces[i]->halfedge() = related_halfedges[i * 2];
         related_halfedges[i * 2]->face() = related_faces[i];
         adjacent_halfedges[i]->face() = related_faces[i];
         related_halfedges[((i + 1) % 4) * 2 + 1]->face() = related_faces[i];
         
         adjacent_vertices[i]->halfedge() = adjacent_halfedges[i]; //related_halfedges[i * 2 + 1];

         related_edges[i]->halfedge() = related_halfedges[i * 2];
      }




      return midpoint;
	 }

   VertexIter HalfedgeMesh::collapseEdge( EdgeIter e )
   {
      // TODO This method should collapse the given edge and return an iterator to the new vertex created by the collapse.

			return VertexIter();
   }

   EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
   {
      // TODO This method should flip the given edge and return an iterator to the flipped edge.

      // get related halfedges
      HalfedgeIter h = e0->halfedge();
      if (h->isBoundary()) {
         return e0;
      }
      HalfedgeIter first_h = h->next();
      HalfedgeIter second_h = first_h->next();

      HalfedgeIter h_twin = h->twin();
      HalfedgeIter first_h_twin = h_twin->next();
      HalfedgeIter second_h_twin = first_h_twin->next();


      // change vertex and vertex degree
         // vertexDegree[h->vertex()]--;
         // vertexDegree[h_twin->vertex()]--;

      h->vertex()->halfedge() = first_h_twin;
      h->face()->halfedge() = first_h_twin;
      h_twin->vertex()->halfedge() = first_h;
      h_twin->face()->halfedge() = first_h;

      h->vertex() = first_h->next()->vertex();
      h_twin->vertex() = first_h_twin->next()->vertex();
         // vertexDegree[h->vertex()]++;
         // vertexDegree[h_twin->vertex()]++;


      // change next
      h->next() = second_h_twin;
      second_h_twin->next() = first_h;
      first_h->next() = h;

      h_twin->next() = second_h;
      second_h->next() = first_h_twin;
      first_h_twin->next() = h_twin;

      printf("!!!!\n");
      printf("%p\n", &e0->halfedge()->vertex());
      return EdgeIter();

   }

   void MeshResampler::upsample( HalfedgeMesh& mesh )
   // This routine should increase the number of triangles in the mesh using Loop subdivision.
   {
      // Each vertex and edge of the original surface can be associated with a vertex in the new (subdivided) surface.
      // Therefore, our strategy for computing the subdivided vertex locations is to *first* compute the new positions
      // using the connectity of the original (coarse) mesh; navigating this mesh will be much easier than navigating
      // the new subdivided (fine) mesh, which has more elements to traverse.  We will then assign vertex positions in
      // the new mesh based on the values we computed for the original mesh.


      // TODO Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
      // TODO and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
      // TODO a vertex of the original mesh.


      // TODO Next, compute the updated vertex positions associated with edges, and store it in Edge::newPosition.


      // TODO Next, we're going to split every edge in the mesh, in any order.  For future
      // TODO reference, we're also going to store some information about which subdivided
      // TODO edges come from splitting an edge in the original mesh, and which edges are new,
      // TODO by setting the flat Edge::isNew.  Note that in this loop, we only want to iterate
      // TODO over edges of the original mesh---otherwise, we'll end up splitting edges that we
      // TODO just split (and the loop will never end!)


      // TODO Now flip any new edge that connects an old and new vertex.


      // TODO Finally, copy the new vertex positions into final Vertex::position.
   }

   // Given an edge, the constructor for EdgeRecord finds the
   // optimal point associated with the edge's current quadric,
   // and assigns this edge a cost based on how much quadric
   // error is observed at this optimal point.
   EdgeRecord::EdgeRecord( EdgeIter& _edge )
   : edge( _edge )
   {
      // TODO Compute the combined quadric from the edge endpoints.


      // TODO Build the 3x3 linear system whose solution minimizes
      // the quadric error associated with these two endpoints.


      // TODO Use this system to solve for the optimal position, and
      // TODO store it in EdgeRecord::optimalPoint.


      // TODO Also store the cost associated with collapsing this edge
      // TODO in EdgeRecord::Cost.

   }

   void MeshResampler::downsample( HalfedgeMesh& mesh )
   {
      // TODO Compute initial quadrics for each face by simply writing the plane
      // equation for the face in homogeneous coordinates.  These quadrics should
      // be stored in Face::quadric


      // TODO Compute an initial quadric for each vertex as the sum of the quadrics
      // associated with the incident faces, storing it in Vertex::quadric


      // TODO Build a priority queue of edges according to their quadric error cost,
      // TODO i.e., by building an EdgeRecord for each edge and sticking it in the queue.


      // TODO Until we reach the target edge budget, collapse the best edge.  Remember
      // TODO to remove from the queue any edge that touches the collapsing edge BEFORE
      // TODO it gets collapsed, and add back into the queue any edge touching the collapsed
      // TODO vertex AFTER it's been collapsed.  Also remember to assign a quadric to the
      // TODO collapsed vertex, and to pop the collapsed edge off the top of the queue.
   }

   void Vertex::computeCentroid( void )
   {
      // TODO Compute the average position of all neighbors of this vertex, and
      // TODO store it in Vertex::centroid.  This value will be used for resampling.
   }

   Vector3D Vertex::normal( void ) const
   // TODO Returns an approximate unit normal at this vertex, computed by
   // TODO taking the area-weighted average of the normals of neighboring
   // TODO triangles, then normalizing.
   {
      // TODO Compute and return the area-weighted unit normal.
			return Vector3D();
	 }

   void MeshResampler::resample( HalfedgeMesh& mesh )
   {
      // TODO Compute the mean edge length.


      // TODO Repeat the four main steps for 5 or 6 iterations


      // TODO Split edges much longer than the target length (being careful about how the loop is written!)


      // TODO Collapse edges much shorter than the target length.  Here we need to be EXTRA careful about
      // TODO advancing the loop, because many edges may have been destroyed by a collapse (which ones?)

      //
      // TODO Now flip each edge if it improves vertex degree


      // TODO Finally, apply some tangential smoothing to the vertex positions
   }
}
