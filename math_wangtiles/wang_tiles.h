
#ifndef DEF_WANGTILES
	#define DEF_WANGTILES

	#include "vec.h"

	struct Tile
	{
		int n, e, s, w;
		int numSubtiles, numSubdivs, numPoints, numSubPoints;
		int ** subdivs;
		Vector2DF *points;
		Vector2DF *subPoints;
	};

	class WangTiles {
	public:
		WangTiles ();
		
		bool	LoadTileSet ( const char* fileName );
		void	SetDensityFunc ( float* density, int xres, int yres );
		void	SetMaxPoints (int m );

		int		RecurseTileImage ( Vector2DF cmin, Vector2DF cmax, float zm, float ts );
		void	RecurseTileImage (Tile & t, float x, float y, int level);

		int numPnts ()				{ return mNumPnts; }
		Vector3DF getPnt ( int n )	{ return mPoints[n]; }

	private:
		
		float		mZoom;
		Vector2DF	mClipMin, mClipMax;

		float*		mDensity;				// input density function
		int			mXRes, mYRes;
	
		int			mNumPnts, mMaxPnts;		// output points 
		Vector3DF*	mPoints;
		Vector3DF*  mCurrPnt;

		Tile*		mTiles;					// wang tile data
		int			numTiles;
		int			numSubtiles;
		int			numSubdivs;

		float		toneScale;
	};

#endif


