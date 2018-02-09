#pragma once
#include <TriMesh.h>
#include <TriMesh_algo.h>
#include <Volume.h> //Tao's volume
#include <spaceinfo.h>
#include <geomalgo.h>

#include <functional>
#include <random>
#include <limits>
#include <cmath>

typedef trimesh::TriMesh TriMesh;
typedef trimesh::TriMesh::Face TriFace;
using trimesh::point;
typedef trimesh::Vec<3, double> dvec3;
typedef trimesh::Vec<2, double> dvec2;
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
#define _USE_MATH_DEFINES
	// epsilon
	double eps = 0.000000001;
	// compute determinant of 3x3 (upper-left-sub-)matrix (input must be at least 3x3)
	double det3x3( const xform& _m );

	// the function: a triangle -> (intensity,weight)
	// used to derive an intensity volume given a mesh
	using std::function;
	typedef function<vec2( TriMesh*, int, const point& )> intenFunc;
	struct SimpleInten
	{
		vec2 operator () (TriMesh* _msh, int _fi, const point& _v)
		{
			auto c = _msh->centroid( _fi );
			vec3 v2c = c - _v; 
			trimesh::normalize( v2c );
			vec3 f_nml = _msh->trinorm( _fi );
			trimesh::normalize( f_nml );
			float scalar = f_nml.dot( v2c );
			float w = 1.0 / std::max( (double)trimesh::dist2( c, _v ), eps );
			return vec2( scalar*w, w );
		}
	};
	struct MVCInten
	{
		inline int i_m_1( int _i ) { return ( _i - 1 + 3 ) % 3; }
		inline int i_p_1( int _i ) { return ( _i + 1 ) % 3; }
		vec2 operator() ( TriMesh* _msh, int _fi, const point& _x )
		{
			double cur_eps = _msh->bbox.size().max() * eps;
			auto f = _msh->faces[ _fi ];
			double w = 0.0, scalar = 0.0;
			vec3 u[ 3 ];
			double d[ 3 ];
			for ( auto i = 0; i < 3; ++i )
			{
				auto pi = _msh->vertices[ f[i] ];
				d[i] = trimesh::dist( pi, _x );
				if ( d[i] < cur_eps )
				{
					w = numeric_limits<float>::max();
					scalar = 1.0f;
					goto MVCIntenRET;
				}
				u[ i ] = ( pi - _x ) / (float)d[ i ];
			}
			double theta[ 3 ];
			for ( auto i = 0; i < 3; ++i )
			{
				auto li = trimesh::len( u[ i_p_1( i ) ] - u[ i_m_1( i ) ] );
				theta[i] = 2 * asin( li / 2.0 );
			}
			auto h = ( theta[ 0 ] + theta[ 1 ] + theta[ 2 ] ) / 2.0f;
			if ( util::is_equal( M_PI, (double)h, (double)eps ) )
			{
				// x lies on this face, use 2d barycentric
				for ( auto i = 0; i < 3; ++i )
				{
					auto wi = std::sin( theta[ i ] ) * d[ i_m_1( i ) ] * d[ i_p_1( i ) ];
					scalar += wi;
				}
				w = numeric_limits<float>::max();
				goto MVCIntenRET;
			}
			{
				double c[ 3 ], s[ 3 ];
				xform orient_mat(
					u[ 0 ][ 0 ], u[ 0 ][ 1 ], u[ 0 ][ 2 ], 0.f,
					u[ 1 ][ 0 ], u[ 1 ][ 1 ], u[ 1 ][ 2 ], 0.f,
					u[ 2 ][ 0 ], u[ 2 ][ 1 ], u[ 2 ][ 2 ], 0.f,
					0.f, 0.f, 0.f, 1.f );
				double sign_orient = det3x3( orient_mat ) > 0.0 ? 1 : -1;
				for ( auto i = 0; i < 3; ++i )
				{
					c[ i ] = ( 2 * std::sin( h ) * std::sin( h - theta[ i ] ) ) /
						( std::sin( theta[ i_p_1( i ) ] ) * std::sin( theta[ i_m_1( i ) ] ) )
						- 1.0;
					s[ i ] = sign_orient * std::sqrt( std::max( eps, 1.0 - c[ i ] * c[ i ] ) );
					if ( util::is_equal( s[ i ], 0.0, cur_eps ) )
					{
						// ignore this face since x on same plane but is outside
						scalar = 0.0;
						w = 0.0;
						goto MVCIntenRET;
					}
				}
				for ( auto i = 0; i < 3; ++i )
				{
					double wi = ( theta[ i ] - c[ i_p_1( i ) ] * theta[ i_m_1( i ) ] - c[ i_m_1( i ) ] * theta[ i_p_1( i ) ] ) /
						2.0 * ( d[ i ] * std::sin( theta[ i_p_1( i ) ] )*s[ i_m_1( i ) ] );
					scalar += wi;
					w += wi;
				}
			}
		MVCIntenRET:
			return vec2( scalar, w );
		}
	};
	struct SDFInten
	{
		vec2 operator() ( TriMesh* _msh, int _fi, const point& _x )
		{
			double d = 0.0;
			double good_face = 1; // good by default
			auto f = _msh->faces[ _fi ];

			// dist from x to the plane of face fi
			vec3 nrml = _msh->trinorm( _fi );
			double nrml_l = trimesh::len( nrml );
			if ( nrml_l < eps )
			{
				d = 0.0;
			}
			else
			{
				nrml /= nrml_l;
				vec3 x2v = _msh->vertices[ f[ 0 ] ] - _x;
				d = x2v.dot( nrml );
			}

			/* does x project to inside of the face? 
			** compute normal of sub-tri: <proj, v_i, v_i+1>,
			** then check consistency with the input tri normal
			*/
			auto proj_x = _x + float( d ) * nrml;
			vec3 x2v[ 3 ];
			for ( auto i = 0; i < 3; ++i )
				x2v[ i ] = _msh->vertices[ f[ i ] ] - _x;
			// compute sub-tri normal
			vec3 nrml_subtri;
			for ( auto i = 0; i < 3; ++i )
			{
				nrml_subtri = x2v[ i ].cross( x2v[ ( i + 1 ) % 3 ] );
				if ( nrml_subtri.dot( nrml ) < 0.0 ) // inconsistent?
				{
					good_face = -1;
					break;
				}
			}
			return vec2( std::abs(d), good_face );
		}
	};
	/*
	** implementation of common defs
	*/
	double det3x3( const xform& _m )
	{
		return 
			_m( 0, 0 ) * ( _m( 1, 1 )*_m( 2, 2 ) - _m( 2, 1 )*_m( 1, 2 ) ) -
			_m( 1, 0 ) * ( _m( 0, 1 )*_m( 2, 2 ) - _m( 2, 1 )*_m( 0, 2 ) ) +
			_m( 2, 0 ) * ( _m( 0, 1 )*_m( 1, 2 ) - _m( 1, 1 )*_m( 0, 2 ) );
	}

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
		static void mkSDFVol( TriMesh* _msh, const trimesh::xform& _vox2mesh, Volume* _vol, const intenFunc& _ifun );
		/*
		** operations on volume
		*/
		static void addNoise( Volume* _vol, double _g_v );
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
		// flip normals & face-orientation if closest point's neighbor faces' normals are not outward
		_msh->need_adjacentfaces();
		const auto& adj_faces = _msh->adjacentfaces[ near_vi ];
		auto near_v = _msh->vertices[ near_vi ];
		bool good_orient = false;
		vec3 nrml = _msh->normals[ near_vi ];
		if ( ( outside_p - near_v ).dot( nrml ) > 0.0f )
		{
			good_orient = true;
		}
		/*for ( auto fi : adj_faces )
		{
			vec3 nrml = _msh->trinorm( fi );
			if ( ( outside_p - near_v ).dot( nrml ) > 0.0f )
			{
				good_orient = true;
				break;
			}
		}*/
		if ( !good_orient )
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
							/*if ( ii == 0 && jj == 0 && kk == 2 && fi == 472 )
								int stop = 1;*/
							auto f = faces[ fi ];
							auto i_res = _ifun( _msh, fi, v );
							if ( i_res[ 1 ] == numeric_limits<float>::max() )
							{
								intensity = i_res[ 0 ];
								w = 1.0;
								break;
							}
							intensity += i_res[0];
							w += i_res[ 1 ];
						}
						//intensity /= std::max( w, eps );
						if ( intensity != intensity )
							printf( "warning: intensity(%d,%d,%d)=%f\n", ii, jj, kk, intensity );
						_vol->setDataAt( ii, jj, kk, intensity );
					}
				}
			}
	}
	void Generator::mkSDFVol( TriMesh* _msh, const trimesh::xform& _vox2mesh, Volume* _vol, const intenFunc& _ifun )
	{
		const auto& faces = _msh->faces;
		const auto& vts = _msh->vertices;
		point vox_c, v_smpl;
		double min_val = numeric_limits<double>::max();
		double max_val = numeric_limits<double>::min();
		for ( int i = 0; i < _vol->getSizeX(); ++i )
		{
			for ( int j = 0; j < _vol->getSizeY(); ++j )
			{
				for ( int k = 0; k < _vol->getSizeZ(); ++k )
				{
					vox_c = point( i, j, k ) + point( 0.5f );
					SpaceConverter::fromVoxToModel( vox_c, v_smpl, _vox2mesh );
					
					// dbg
					double min_d = numeric_limits<double>::max();
					int min_vi = 0;
					for ( auto ii = 0; ii < vts.size(); ++ii )
					{
						auto v = vts[ ii ];
						auto d = trimesh::dist( v, v_smpl );
						if ( d < min_d )
						{
							min_d = d;
							min_vi = ii;
						}
					}
					double sign = _msh->normals[ min_vi ].dot( vts[ min_vi ] - v_smpl ) > 0.0 ? 1 : -1;
					double d = sign * min_d;

					/* now compute signed distance to the closest element (face / vert) */
					// first, find the distance to the closest mesh face
					//double min_d = numeric_limits<double>::max();
					//int fi_cand = -1;
					//for ( auto fi = 0; fi < vts.size(); ++fi )
					//{
					//	auto res = _ifun( _msh, fi, v_smpl );
					//	auto d = res[ 0 ];
					//	auto good_face = res[ 1 ];
					//	if ( good_face > 0.0f && d < min_d )
					//	{
					//		fi_cand = fi;
					//		min_d = d;
					//	}
					//}
					//auto f = faces[ fi_cand ];
					//// second, iterate thru each vert to see if even smaller distance can be found
					//int min_vi = -1;
					//for ( auto vi = 0; vi < vts.size(); ++vi )
					//{
					//	auto v = vts[ vi ];
					//	auto d = trimesh::dist( v, v_smpl );
					//	if ( d < min_d )
					//	{
					//		min_d = d;
					//		min_vi = vi;
					//	}
					//}
					//int vi = min_vi >= 0 ? min_vi : faces[ fi_cand ][ 0 ];
					//// the sign of the distance
					//double sign = _msh->normals[ vi ].dot( vts[ vi ] - v_smpl ) > 0.0 ? 1 : -1;
					//// now we have the signed-distance
					//double d = sign*min_d;
					_vol->setDataAt( i, j, k, d );
					min_val = std::min( min_val, d );
					max_val = std::max( max_val, d );
				}
			}
		}
		printf( "SDF vol range [%f, %f] \n", min_val, max_val );
	}
	void Generator::addNoise( Volume* _vol, double _g_v )
	{
		// use constant seed for reproducible random sequence
		unsigned seed = 1000; //std::chrono::system_clock::now().time_since_epoch().count();
		std::default_random_engine generator( seed );
		std::normal_distribution<double> guassian( 0.0, _g_v );
		std::uniform_real_distribution<double> uniform( 0.0, 1.0 );
		
		double noise_min = numeric_limits<double>::max();
		double noise_max = numeric_limits<double>::min();
		for ( auto i = 0; i < _vol->getSizeX(); ++i )
		{
			for ( auto j = 0; j < _vol->getSizeY(); ++j )
			{
				for ( auto k = 0; k < _vol->getSizeZ(); ++k )
				{
					double old_v = _vol->getDataAt( i, j, k );
					double noise = guassian( generator );
					// Russian roulette to randomly add/discard noise
					bool discard = uniform( generator ) > 0.5;
					if ( discard )
						continue;
					noise_min = std::min( noise_min, noise );
					noise_max = std::max( noise_max, noise );
					_vol->setDataAt( i, j, k, old_v + noise );
				}
			}
		}
		printf( "added noise from range [%f, %f] \n", noise_min, noise_max );
	}
	xform Generator::computeSpaceTransform( const box3& _bbox1, const box3& _bbox2 )
	{
		auto c1 = _bbox1.center();
		auto c2 = _bbox2.center();
		// find max len of box2
		int max_axis = 0;
		double s2 = 0.0;
		for ( auto i = 0; i < 3; ++i )
		{
			if ( s2 < _bbox2.size()[ i ] )
			{
				max_axis = i;
				s2 = _bbox2.size()[ i ];
			}
		}
		s2 *= 1.1; // shrink mesh a little bit
		double r = s2 / _bbox1.size()[max_axis];
		return xform::trans( c2 ) * xform::scale( r ) * xform::trans( -c1 );
	}
}