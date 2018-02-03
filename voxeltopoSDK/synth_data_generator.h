#pragma once
#include <TriMesh.h>
#include <Volume.h> //Tao's volume
#include <functional>

typedef trimesh::TriMesh TriMesh;
typedef trimesh::TriMesh::Face TriFace;
using trimesh::point;
using trimesh::vec3;
using trimesh::ivec3;
using trimesh::ivec2;
using trimesh::box3;

namespace synthdata
{
	// the function: a triangle -> (intensity,weight)
	// used to derive an intensity volume given a mesh
	using std::function;
	typedef function<float( vec3, vec3, vec3, vec3, vec3, vec3 )> intenFunc;
	class Generator
	{
	private:
		Generator() {}
		Generator( const Generator& ) {}
		Generator& operator = ( const Generator& ) {}
	protected:
		static point mk_meshspace_sample( const Volume* _vol, int _ii, int _jj, int _kk, 
			const point& _o, const box3& _bbox );
	public:
		/*
		** make a dense volume from given triangle mesh
		*/
		// make an intensity volume
		static void mkIntensityVol( const TriMesh* _msh, Volume* _vol, const intenFunc& _ifun );
		// make a signed-distance-function volume
		static void mkSDFVol( const TriMesh* _msh, Volume* _vol );
	};

	void Generator::mkIntensityVol( const TriMesh* _msh, Volume* _vol, const intenFunc& _ifun )
	{
		const auto &faces = _msh->faces;
		for ( auto fi = 0; fi < faces.size(); ++fi )
		{
			auto f = faces[ fi ];
			for ( auto ii = 0; ii < _vol->getSizeX(); ++ii )
			{
				for ( auto jj = 0; jj < _vol->getSizeY(); ++jj )
				{
					for ( auto kk = 0; kk < _vol->getSizeX(); ++kk )
					{
						point v = mk_meshspace_sample( ii, jj, kk, o, bbox );
					}
				}
			}
		}
	}
	void Generator::mkSDFVol( const TriMesh* _msh, Volume* _vol )
	{

	}
}