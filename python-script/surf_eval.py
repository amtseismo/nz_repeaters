#!/usr/bin/env python

## @file surf_eval.py

## @brief Python application to evaluate bivariate b-spline at a set of points,
## given the surface parameters and point coordinates. Coordinate
## transformations are also performed.

import math
import numpy
import scipy.interpolate
import shapely.geometry
import sys
# import pdb
from pyproj import Proj
from pyproj import transform

from pyre.applications.Script import Script as Application

class SurfEval(Application):
  """
  Python application to evaluate bivariate b-spline at a set of points,
  given the surface parameters and point coordinates. Coordinate
  transformations are also performed.
  """

  ## \b Properties
  ## @li \b surf_param_file Input file with surface definition parameters.
  ## @li \b sample_method Method of sampling.
  ## @li \b sample_x_min Min X value for rectilinear sampling.
  ## @li \b sample_y_min Min Y value for rectilinear sampling.
  ## @li \b sample_x_max Max X value for rectilinear sampling.
  ## @li \b sample_y_max Max Y value for rectilinear sampling.
  ## @li \b sample_num_x Number of X points for rectilinear sampling.
  ## @li \b sample_num_y Number of Y points for rectilinear sampling.
  ## @li \b sample_coord_file Sample filename for from_file sampling.
  ## @li \b coord_convert Perform coordinate conversion?
  ## @li \b valid_region Polygon file defining valid sampling region.
  ## @li \b out_of_range_option Whether to exclude points or use a fill value.
  ## @li \b out_of_range_fill_val Fill value used for points that are out of range.
  ## @li \b output_txt Whether to create tab-delimited text output file.
  ## @li \b output_txt_file Output file with surface coordinates.
  ## @li \b output_vtk Whether to create output VTK file.
  ## @li \b output_vtk_file Output VTK file with surface coordinates.
  ## @li \b sample_coordsys Proj parameters defining sample coordinate system.
  ## @li \b z_scale Scaling factor for interface elevation (relative to m).
  ## @li \b compute_fit Compute fit to z-values at sampling points.
  ## @li \b compute_normal Compute normal vector at sampling points.

  import pyre.inventory
  
  surfParamFile = pyre.inventory.str("surf_param_file",
                                     default="surf_params.npz")
  surfParamFile.meta['tip'] = "Input file with surface definition parameters."
  
  sampleMethod = pyre.inventory.str(
    "sample_method", default="rectilinear",
    validator=pyre.inventory.choice(["rectilinear","from_file","from_array"]))
  sampleMethod.meta['tip'] = "Method of specifying sample points."
  
  sampleXMin = pyre.inventory.float("sample_x_min", default=-50000.0)
  sampleXMin.meta['tip'] = "Min X value for rectilinear sampling."
  
  sampleXMax = pyre.inventory.float("sample_x_max", default=50000.0)
  sampleXMax.meta['tip'] = "Max X value for rectilinear sampling."
  
  sampleNumX = pyre.inventory.int("sample_num_x", default=11)
  sampleNumX.meta['tip'] = "Number of X points for rectilinear sampling."
  
  sampleYMin = pyre.inventory.float("sample_y_min", default=-50000.0)
  sampleYMin.meta['tip'] = "Min Y value for rectilinear sampling."
  
  sampleYMax = pyre.inventory.float("sample_y_max", default=50000.0)
  sampleYMax.meta['tip'] = "Max Y value for rectilinear sampling."
  
  sampleNumY = pyre.inventory.int("sample_num_y", default=11)
  sampleNumY.meta['tip'] = "Number of Y points for rectilinear sampling."
  
  sampleCoordFile = pyre.inventory.str("sample_coord_file",
                                       default="sample_coords.txt")
  sampleCoordFile.meta['tip'] = "Filename with coords for from_file sampling."
  
  coordConvert = pyre.inventory.bool("coord_convert", default=True)
  coordConvert.meta['tip'] = "Whether to convert local coordinates to sample CS."
  
  validRegion = pyre.inventory.str("valid_region", default="valid_region.txt")
  validRegion.meta['tip'] = "Polygon file defining valid sampling region."
  
  outOfRangeOption = pyre.inventory.str("out_of_range_option",
                                        default="exclude",
                                        validator=pyre.inventory.choice(
      ["exclude","fill_with_fill_val"]))
  outOfRangeOption.meta['tip'] = "Whether to exclude out of range values or use fill value."

  outOfRangeFillVal = pyre.inventory.str("out_of_range_fill_val",
                                         default="nan")
  outOfRangeFillVal.meta['tip'] = "Fill value for excluded values."
  
  outputAsArray = pyre.inventory.bool("output_as_array", default=False)
  outputAsArray.meta['tip'] = "Whether to create an output array attribute for this class, that can serve as an interface to other programs."

  outputTxt = pyre.inventory.bool("output_txt", default=True)
  outputTxt.meta['tip'] = "Whether to create column-delimited output file."
  
  outputTxtFile = pyre.inventory.str("output_txt_file",
                                     default="surface_points.txt")
  outputTxtFile.meta['tip'] = "Output file with surface coordinates."
  
  outputVtk = pyre.inventory.bool("output_vtk", default=False)
  outputVtk.meta['tip'] = "Whether to create column-delimited output file."
  
  outputVtkFile = pyre.inventory.str("output_vtk_file",
                                     default="surface_points.vtk")
  outputVtkFile.meta['tip'] = "Output VTK file with surface coordinates."
  
  sampleCoordsys = pyre.inventory.str(
    "sample_coordsys",
    default="+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000 +y_0=10000000 +a=6378137.0 +rf=298.257222101 +towgs84=0.0,0.0,0.0")
  sampleCoordsys.meta['tip'] = "Proj parameters for sample coordinate system."
  
  zScale = pyre.inventory.float("z_scale", default=1.0e-3)
  zScale.meta['tip'] = "Scaling factor for interface elevation (relative to m)."
  
  computeFit = pyre.inventory.bool("compute_fit", default=False)
  computeFit.meta['tip'] = "Compute fit to z-value at sampling points."
  
  computeNormal = pyre.inventory.bool("compute_normal", default=False)
  computeNormal.meta['tip'] = "Compute surface normal at sampling points."

  beVerbose = pyre.inventory.bool("give_debugging_output", default=True)
  beVerbose.meta['tip'] = "Print information about the program progress and diagnostic information"


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="surf_eval"):
    Application.__init__(self, name)
    self.numSamples = 0
    self.numCropped = 0
    self.numValidSamples = 0
    self.sampleCoords = None
    self.sampleCoordsFit = None
    self.sampleCoordsOut = None
    self.normalVec = None
    self.projSample = None
    self.projFit = None
    self.misfit = None
    self.validSample = []
    self.surfParams = []
    self.fitCoordsys = \
      "+proj=tmerc +lon_0=175.45 +lat_0=-40.825 +ellps=WGS84 +datum=WGS84 +k=0.9996 +towgs84=0.0,0.0,0.0"
    self.fitRot = math.radians(45.0)

    return


  def main(self):
    # pdb.set_trace()
    if self.beVerbose: print ""
    self.fillVal = float(self.outOfRangeFillVal)
    self._getSampleCoords()
    self._readSurfParams()
    self._readPolygons()
    self._getSurfCoords()
    self._writeSurfCoords()

    return


  # PRIVATE METHODS /////////////////////////////////////////////////////

  def _getCoordinatesFromNumpyArray(self, coordinatesAsArray):
      self.coordinatesFromArray = coordinatesAsArray
      
      
  def _getSampleCoords(self):
    """
    Function to get sample coordinates.
    """

    if (self.sampleMethod == "from_file"):
      sampleCoords = numpy.loadtxt(self.sampleCoordFile, dtype=numpy.float64)
      if self.beVerbose: print "Sample coordinate file:   %s" % self.sampleCoordFile
      sys.stdout.flush()
      xCoords = sampleCoords[:,0]
      yCoords = sampleCoords[:,1]
      zCoords = numpy.zeros_like(xCoords)
      if (sampleCoords.shape[1] > 2):
        zCoords = sampleCoords[:,2]
        
      self.sampleCoords = numpy.column_stack((xCoords, yCoords, zCoords))
    elif (self.sampleMethod == "from_array"):
      sampleCoords = self.coordinatesFromArray
      if self.beVerbose: print "Taking Coordinates from provided array"
      sys.stdout.flush()
      xCoords = sampleCoords[:,0]
      yCoords = sampleCoords[:,1]
      zCoords = numpy.zeros_like(xCoords)
      if (sampleCoords.shape[1] > 2):
        zCoords = sampleCoords[:,2]
        
      self.sampleCoords = numpy.column_stack((xCoords, yCoords, zCoords))
      
    else:
      if self.beVerbose: print "Generating rectilinear coordinates"
      xCoords = numpy.linspace(self.sampleXMin, self.sampleXMax,
                               num=self.sampleNumX)
      yCoords = numpy.linspace(self.sampleYMin, self.sampleYMax,
                               num=self.sampleNumY)
      xGrid, yGrid = numpy.meshgrid(xCoords, yCoords)
      zGrid = numpy.zeros_like(xGrid.flatten())
      self.sampleCoords = numpy.column_stack((xGrid.flatten(), yGrid.flatten(),
                                              zGrid))

    self.numSamples = self.sampleCoords.shape[0]

    # Define projections.
    self.projSample = Proj(self.sampleCoordsys)
    self.projFit = Proj(self.fitCoordsys)

    # Convert sample coords to fit coordinate system.
    if (self.coordConvert):
      self.sampleCoordsFit = self._sampleToFit(self.sampleCoords)
    else:
      self.sampleCoordsFit = self.sampleCoords.copy()

    return
    

  def _sampleToFit(self, sampleCoords):
    """
    Function to convert sample coordinates to fit coordinate system.
    """
    
    # Convert to local Transverse Mercator coordinates used for surface fitting.
    (xFit, yFit, zFit) = None, None, None
    xOrig = sampleCoords[:,0]
    yOrig = sampleCoords[:,1]
    zOrig = sampleCoords[:,2]
    (xFit, yFit, zFit) = transform(self.projSample, self.projFit,
                                   xOrig, yOrig, zOrig)

    # Create transpose of coordinate array and rotation matrix.
    coordsRot = numpy.column_stack((xFit, yFit, zFit)).transpose()
    cosRot = math.cos(self.fitRot)
    sinRot = math.sin(self.fitRot)
    rot = numpy.array(([[cosRot, -sinRot, 0.0],
                        [sinRot,  cosRot, 0.0],
                        [0.0,        0.0, 1.0]]), dtype=numpy.float64)

    # Rotate to surface fit coordinate system.
    sampleCoordsFit = numpy.dot(rot, coordsRot).transpose()

    return sampleCoordsFit
    

  def _fitToSample(self, fitCoords):
    """
    Function to convert surface fit coordinates to sample coordinate system.
    """
    
    # Create transpose of coordinate array and rotation matrix.
    coordsTrans = fitCoords.transpose()
    cosRot = math.cos(self.fitRot)
    sinRot = math.sin(self.fitRot)
    rot = numpy.array(([[ cosRot,  sinRot, 0.0],
                        [-sinRot,  cosRot, 0.0],
                        [ 0.0,        0.0, 1.0]]), dtype=numpy.float64)

    # Unrotate surface fit coordinate system.
    coordsRot = numpy.dot(rot, coordsTrans).transpose()

    # Convert unrotated surface fit coords to sample coordinate system.
    (xSample, ySample, zSample) = None, None, None
    xRot = coordsRot[:,0]
    yRot = coordsRot[:,1]
    zRot = coordsRot[:,2]
    (xSample, ySample, zSample) = transform(self.projFit, self.projSample,
                                            xRot, yRot, zRot)
    sampleCoords = numpy.column_stack((xSample, ySample, zSample))

    return sampleCoords


  def _readSurfParams(self):
    """
    Function to read surface definition parameters.
    """

    if self.beVerbose: print "Surface parameter file:   %s" % self.surfParamFile
    data = numpy.load(self.surfParamFile)
    tx = data['tx']
    ty = data['ty']
    c = data['c']
    kx = data['kx']
    ky = data['ky']
    self.surfParams = [tx, ty, c, kx, ky]

    return


  def _readPolygons(self):
    """
    Function to read polygon definitions and determine which points lie within
    each polygon.
    """

    # Read polygon and create shapely object.
    validPoints = numpy.loadtxt(self.validRegion)
    validRegion = shapely.geometry.Polygon(validPoints)

    # Loop over specified sample points remove those outside valid region.

    self.numCropped = 0
    for pointNum in range(self.numSamples):
      xVal = self.sampleCoordsFit[pointNum,0]
      yVal = self.sampleCoordsFit[pointNum,1]
      point = shapely.geometry.Point(xVal, yVal)
      if (validRegion.intersects(point)):
        self.validSample.append(pointNum)
      else:
        self.numCropped += 1

    self.numValidSamples = len(self.validSample)

    if self.beVerbose: print "Number of original sample points:  %d" % self.numSamples
    if self.beVerbose: print "Number of cropped sample points:   %d" % self.numCropped
    if self.beVerbose: print "Number of valid sample points:     %d" % self.numValidSamples
    if self.beVerbose: print 'Valid region file was: ', self.validRegion
    sys.stdout.flush()

    return

    
  def _getSurfCoords(self):
    """
    Function to transform sample coordinates to coordinates used to define the
    surface, get the surface coordinates, and then transform back to the
    sample coordinate system.
    """

    # Get sample points in surface fit coordinates.
    xFit = self.sampleCoordsFit[self.validSample,0]
    yFit = self.sampleCoordsFit[self.validSample,1]
    zFit = numpy.zeros_like(yFit)

    # Set up info for surface normals.
    xDeriv = numpy.zeros_like(zFit)
    yDeriv = numpy.zeros_like(zFit)
    zDeriv = numpy.ones_like(zFit)

    # Evaluate spline at selected points.
    for pointNum in range(self.numValidSamples):
      x = xFit[pointNum]
      y = yFit[pointNum]
      zFit[pointNum] = scipy.interpolate.bisplev(x, y, self.surfParams)
      if (self.computeNormal):
        xDeriv[pointNum] = scipy.interpolate.bisplev(x, y, self.surfParams,
                                                     dx=1, dy=0)
        yDeriv[pointNum] = scipy.interpolate.bisplev(x, y, self.surfParams,
                                                     dx=0, dy=1)

    # Transform back to sample coordinate system.
    sampleCoordsFit = numpy.column_stack((xFit, yFit, zFit))
    if (self.coordConvert):
      self.sampleCoordsOut = self._fitToSample(sampleCoordsFit)
    else:
      self.sampleCoordsOut = sampleCoordsFit.copy()

    # Compute normals in sample coordinate system, if requested.
    if (self.computeNormal):
      normalVecFit = numpy.column_stack((-xDeriv, -yDeriv, zDeriv))
      vecMag = numpy.apply_along_axis(numpy.linalg.norm, 1,
                               normalVecFit).reshape(-1,1)
      normalVecFit /= vecMag
      sampleCoordsOffFit = sampleCoordsFit + 1000.0 * normalVecFit
      if (self.coordConvert):
        sampleCoordsOff = self._fitToSample(sampleCoordsOffFit)
      else:
        sampleCoordsOff = sampleCoordsOffFit.copy()
      coordDiff = sampleCoordsOff - self.sampleCoordsOut
      diffMag = numpy.apply_along_axis(numpy.linalg.norm, 1,
                                       coordDiff).reshape(-1,1)
      self.normalVec = coordDiff/diffMag
                                                       
    # Scale z-values.
    xSample = self.sampleCoordsOut[:,0]
    ySample = self.sampleCoordsOut[:,1]
    zSample = self.sampleCoordsOut[:,2] * self.zScale
    self.sampleCoordsOut[:,2] = zSample

    # Output some info about computed coordinates.
    maxInd = numpy.argmax(zSample)
    minInd = numpy.argmin(zSample)
    xMax = xSample[maxInd]
    xMin = xSample[minInd]
    yMax = ySample[maxInd]
    yMin = ySample[minInd]
    zMax = zSample[maxInd]
    zMin = zSample[minInd]
    posVals = numpy.nonzero(zSample > 0.0)
    numPos = posVals[0].shape[0]
    if self.beVerbose: print "Maximum elevation value:        %g, %g, %g" % (xMax, yMax, zMax)
    if self.beVerbose: print "Minimum elevation value:        %g, %g, %g" % (xMin, yMin, zMin)
    if self.beVerbose: print "Number of positive elevations:  %d" % numPos
    sys.stdout.flush()

    # Compute misfit, if requested.
    if (self.computeFit):
      self.misfit = (self.sampleCoords[self.validSample,2] - zSample).reshape((
        self.numValidSamples, 1))
      norm = numpy.linalg.norm(self.misfit)
      normAvg = norm/float(self.numValidSamples)
      minInd = numpy.argmin(numpy.abs(self.misfit))
      maxInd = numpy.argmax(numpy.abs(self.misfit))
      minMisfit = abs(self.misfit[minInd])
      minX = xSample[minInd]
      minY = ySample[minInd]
      maxMisfit = abs(self.misfit[maxInd])
      maxX = xSample[maxInd]
      maxY = ySample[maxInd]
      if self.beVerbose: print "Norm of elevation misfit:           %g" % norm
      if self.beVerbose: print "Average norm of elevation misfit:   %g" % normAvg
      if self.beVerbose: print "Minimum absolute misfit:            %g" % minMisfit
      if self.beVerbose: print "  Coordinates:                      %g, %g" % (minX, minY)
      if self.beVerbose: print "Maximum absolute misfit:            %g" % maxMisfit
      if self.beVerbose: print "  Coordinates:                      %g, %g" % (maxX, maxY)
      sys.stdout.flush()

    return

    
  def _writeSurfCoords(self):
    """
    Function to write out text and VTK files with the surface coordinates.
    """
      
    # Set up output arrays, depending on whether we want to exclude points or
    # use fill values.
    if (self.outOfRangeOption == "exclude"):
      coordsOut = self.sampleCoordsOut
      misfitOut = self.misfit
      normalOut = self.normalVec
    else:
      coordsOut = numpy.zeros((self.numSamples, 3))
      coordsOut[:,2] = self.fillVal
      coordsOut[:,0:2] = self.sampleCoords[:,0:2]
      coordsOut[self.validSample,:] = self.sampleCoordsOut
      misfitOut = self.fillVal * numpy.ones((self.numSamples))
      if (self.computeFit):
        misfitOut[self.validSample] = self.misfit
      normalOut = self.fillVal * numpy.ones((self.numSamples, 3))
      if (self.computeNormal):
        normalOut[self.validSample, :] = self.normalVec
      
    numPointsOut = coordsOut.shape[0]

    # Create output array, if requested.
    if (self.outputAsArray):
      self.OutputArray = coordsOut
      if (self.computeFit):
        self.OutputArray = numpy.hstack((self.OutputArray, misfitOut))
      if (self.computeNormal):
        self.OutputArray = numpy.hstack((self.OutputArray, normalOut))
        
    if (self.outputTxt):
      outVals = coordsOut
      head = 'X-Coord\tY-Coord\tZ-Coord'
      if (self.computeFit):
        outVals = numpy.hstack((outVals, misfitOut))
        head += '\tZ-Misfit'
      if (self.computeNormal):
        outVals = numpy.hstack((outVals, normalOut))
        head += '\tX-Normal\tY-Normal\tZ-Normal'
      
        
      head += '\n'
      f = open(self.outputTxtFile, 'w')
      f.write(head)
      numpy.savetxt(f, outVals, delimiter='\t')
      f.close()

    # Output VTK file.
    if (self.outputVtk):
      wPoint = False
      vtkHead = "# vtk DataFile Version 2.0\n" + \
        "Surface coordinates\n" + \
        "ASCII\n" + \
        "DATASET POLYDATA\n" + \
        "POINTS " + repr(numPointsOut) + " double\n"

      v = open(self.outputVtkFile, 'w')
      v.write(vtkHead)
      numpy.savetxt(v, coordsOut)

      if (self.computeFit):
        sHead = "POINT_DATA " + repr(numPointsOut) + "\n" + \
          "SCALARS obs_minus_pred double 1\n" + \
          "LOOKUP_TABLE DEFAULT\n"
        wPoint = True
        v.write(sHead)
        numpy.savetxt(v, misfitOut)

      if (self.computeNormal):
        if (wPoint):
          vHead = "VECTORS normal_vector double\n"
        else:
          vHead = "POINT_DATA " + repr(numPointsOut) + "\n" + \
            "VECTORS normal_vector double\n"
        v.write(vHead)
        numpy.savetxt(v, normalOut)
      
      v.close()

    return

# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = SurfEval()
  app.run()

# End of file
