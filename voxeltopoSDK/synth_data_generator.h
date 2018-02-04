#pragma once
#include <TriMesh.h>
#include <TriMesh_algo.h>
#include <Volume.h> //Tao's volume
#include <spaceinfo.h>
#include <functional>
#include <limits>

typedef trimesh::TriMesh TriMesh;
typedef trimesh::TriMesh::Face TriFace;
using trimesh::point;
using trimesh::vec3;
using trimesh::vec2;
using trimesh::ivec3;
using trimesh::ivec2;
using trimesh::box3;
using trimesh::xform;
using std::numeric_limits;

namespace synthdata
{
	/*
	** common defs
	*/
	// epsilon
	double eps = 0.000000001;

	// the function: a triangle -> (intensity,weight)
	// used to derive an intensity volume given a mesh
	using std::function;
	typedef function<vec2( TriMesh*, int, point )> intenFunc;
	struct SimpleInten
	{
		vec2 operator () (TriMesh* _msh, int _fi, point _v)
		{
			auto c = _msh->centroid( _fi );
			vec3 v2c = c - _v; 
			trimesh::normalize( v2c );
			vec3 f_nml = _msh->trinorm( _fi );
			trimesh::normalize( f_nml );
			float scalar = f_nml.dot( v2c );
			float w = 1.0 / std::max( (double)trimesh::dist2( c, _v ), eps );
			return vec2( scalar, w );
		}
	};
	struct MVCInten
	{
		vec2 operator() ( TriMesh* _msh, int _fi, point _v )
		{
			return vec2( 1.0f );
		}
	};
	class Generator
	{
	private:
		Generator() {}
		Generator( const Generator& ) {}
		Generator& operator = ( const Generator& ) {}
	protected:
		
	public:
		// compute bbox for given volume
		static box3 mkBBoxForVol( Volume* _vol );
		// transform matrix from space 1 (bbox1) -> space 2 (bbox2)
		static xform computeSpaceTransform( const box3& _bbox1, const box3& _bbox2 );
		// make sure mesh normals point outward
		static void enforceNormalsOutward( TriMesh* _msh );
		/*
		** make a dense volume from given triangle mesh
		*/
		// make an intensity volume
		static void mkIntensityVol( TriMesh* _msh, const trimesh::xform& _vox2mesh, Volume* _vol, const intenFunc& _ifun );
		// make a signed-distance-function volume
		static void mkSDFVol( TriMesh* _msh, Volume* _vol );
	};

	box3 Generator::mkBBoxForVol( Volume* _vol )
	{
		box3 bb;
		bb += point(0.0f);
		bb += point( _vol->getSizeX(), _vol->getSizeY(), _vol->getSizeZ() );
		return bb;
	}
	void Generator::enforceNormalsOutward( TriMesh* _msh )
	{
		// orient all normals consistently
		trimesh::orient( _msh );
		_msh->need_normals();
		// find an outside point 
		_msh->need_bbox();
		auto outside_p = _msh->bbox.min - vec3( _msh->bbox.size().max() );
		// find closest point on mesh
		float near_vi = 0;
		float min_d = numeric_limits<float>::max();
		for ( auto vi = 0; vi < _msh->vertices.size(); ++vi )
		{
			point v = _msh->vertices[ vi ];
			auto d = trimesh::dist( v, outside_p );
			if ( min_d > d )
			{
				min_d = d;
				near_vi = vi;
			}
		}
		// flip normals & face-orientation if closest point's normal is not outward
		vec3 nrml = _msh->normals[ near_vi ];
		if ( ( outside_p - _msh->vertices[ near_vi ] ).dot( nrml ) < 0.0f )
		{
			trimesh::faceflip( _msh );
			for ( auto& n : _msh->normals )
				n *= -1;
		}
	}
	void Generator::mkIntensityVol( TriMesh* _msh, const trimesh::xform& _vox2mesh, Volume* _vol, const intenFunc& _ifun )
	{
		const auto &faces = _msh->faces;
		point v;
		point vox_c;
			for ( auto ii = 0; ii < _vol->getSizeX(); ++ii )
			{
				for ( auto jj = 0; jj < _vol->getSizeY(); ++jj )
				{
					for ( auto kk = 0; kk < _vol->getSizeX(); ++kk )
					{
						vox_c = point( ii, jj, kk ) + point( 0.5f );
						SpaceConverter::fromVoxToModel( vox_c, v, _vox2mesh );
						double intensity = 0.0;
						double w = 0.0;
						// accumulate intensity from each face to this voxel sample
						for ( auto fi = 0; fi < faces.size(); ++fi )
						{
							auto f = faces[ fi ];
							auto i_res = _ifun( _msh, fi, v );
							intensity += i_res[0] * i_res[1];
							w += i_res[ 1 ];
						}
						intensity /= std::max( w, eps );
						_vol->setDataAt( ii, jj, kk, intensity );
					}
				}
			}
	}
	void Generator::mkSDFVol( TriMesh* _msh, Volume* _vol )
	{

	}
	xform Generator::computeSpaceTransform( const box3& _bbox1, const box3& _bbox2 )
	{
		auto c1 = trimesh::mix( _bbox1.min, _bbox1.max, 0.5f );
		auto s1 = _bbox1.size();
		auto c2 = trimesh::mix( _bbox2.min, _bbox2.max, 0.5f );
		auto s2 = _bbox2.size();
		s2 *= 1.1; // shrink mesh a little bit
		auto r = s2 / s1;
		return xform::trans( c2 ) * xform::scale( r[ 0 ], r[ 1 ], r[ 2 ] ) * xform::trans( -c1 );
	}
}