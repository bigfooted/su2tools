"""This file can be loaded as a plugin for ParaView >= 5.6

Author: Bane Sullivan <banesulli@gmail.com>
        rewritten by Nijso Beishuizen for su2 file writing
"""

import numpy as np

# This is module to import. It provides VTKPythonAlgorithmBase, the base class
# for all python-based vtkAlgorithm subclasses in VTK and decorators used to
# 'register' the algorithm with ParaView along with information about UI.
from paraview.util.vtkAlgorithm import *
import vtk
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from paraview.simple import *
from vtkmodules.vtkCommonDataModel import vtkMultiBlockDataSet, vtkUnstructuredGrid, vtkPartitionedDataSet

def createModifiedCallback(anobject):
    """Return a function needed for multiple tag selection"""
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified



# Predefined ranges for values of 'Tag' property
TAG_RANGES = {
    'Dots': (100, 110), # (min, max)
    'Wake points': (110, 110),
    'Polygons': (300, 332),
    'Body Panel': (332, 332),
    }


class WriterBase(VTKPythonAlgorithmBase):
    """This is a writer base class to add convenience methods to the
    ``VTKPythonAlgorithmBase`` for writer algorithms and was originally
    implemented in `PVGeo`_ by `Bane Sullivan`_.

    .. _PVGeo: http://pvgeo.org
    .. _Bane Sullivan: http://banesullivan.com

    For more information on what functionality is available, check out the VTK
    Docs for the `vtkAlgorithm`_ and then check out the following blog posts:

    * `vtkPythonAlgorithm is great`_
    * A VTK pipeline primer `(part 1)`_, `(part 2)`_, and `(part 3)`_
    * `ParaView Python Docs`_

    .. _vtkAlgorithm: https://www.vtk.org/doc/nightly/html/classvtkAlgorithm.html
    .. _vtkPythonAlgorithm is great: https://blog.kitware.com/vtkpythonalgorithm-is-great/
    .. _(part 1): https://blog.kitware.com/a-vtk-pipeline-primer-part-1/
    .. _(part 2): https://blog.kitware.com/a-vtk-pipeline-primer-part-2/
    .. _(part 3): https://blog.kitware.com/a-vtk-pipeline-primer-part-3/
    .. _ParaView Python Docs: https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/paraview.util.vtkAlgorithm.html
    """

    def __init__(self, nInputPorts=1, inputType='vtkPolyData', **kwargs):
        VTKPythonAlgorithmBase.__init__(
            self, nInputPorts=nInputPorts, inputType=inputType, nOutputPorts=0
        )
        self.__filename = kwargs.get('filename', None)
        self.__fmt = '%.9e'
        # For composite datasets: not always used
        self.__blockfilenames = None
        self.__composite = False
        self.__delimiter = ','
        self.__stringval = 'blabla'

        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent",
                                         createModifiedCallback(self))
        self._arrayselection.AddArray("one")
        self._arrayselection.AddArray("two")
        self._arrayselection.AddArray("three")
        self._arrayselection.EnableAllArrays()


    def FillInputPortInformation(self, port, info):
        """Allows us to save composite datasets as well.
        NOTE: I only care about ``vtkMultiBlockDataSet``s but you could hack
        this method and ``RequestData`` to handle ``vtkMultiBlockDataSet``s for
        a general use case.
        """
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), self.InputType)
        info.Append(
            self.INPUT_REQUIRED_DATA_TYPE(), 'vtkMultiBlockDataSet'
        )  # 'vtkCompositeDataSet'
        return 1

    def SetFileName(self, filename):
        """Specify the filename for the output.
        This will be appended if saving composite datasets.
        """
        if not isinstance(filename, str):
            raise RuntimeError(
                'File name must be string. Only single file is supported.'
            )
        if self.__filename != filename:
            self.__filename = filename
            self.Modified()

    def get_file_name(self):
        """Get the set filename."""
        return self.__filename

    def _get_array_selection(self):
        print("selecting the array")
        return self._arrayselection

    def Write(self, input_data_object=None):
        """A Python focused convenience method to perform the write out."""
        if input_data_object:
            self.SetInputDataObject(input_data_object)
        self.Modified()
        self.Update()

    def PerformWriteOut(self, input_data_object, filename, object_name):
        """This method must be implemented. This is automatically called by
        ``RequestData`` for single inputs or composite inputs."""
        raise NotImplementedError('PerformWriteOut must be implemented!')

    def apply(self, input_data_object):
        """A convienace method if using these algorithms in a Python environment."""
        self.SetInputDataObject(input_data_object)
        self.Modified()
        self.Update()

    def set_format(self, fmt):
        """Use to set the ASCII format for the writer default is ``'%.9e'``"""
        if self.__fmt != fmt and isinstance(fmt, str):
            self.__fmt = fmt
            self.Modified()

    def get_format(self):
        return self.__fmt

    def set_delimiter(self, deli):
        """The string delimiter to use"""
        if self.__delimiter != deli:
            self.__delimiter = deli
            self.Modified()

    def set_string(self, stringval):
        """The string to use"""
        if self.__stringval != stringval:
            self.__stringval = stringval
            self.Modified()

    def set_array(self, arrayval):
        """The string to use"""
        self._arrayselection = vtkDataArraySelection()
        self._arrayselection.AddObserver("ModifiedEvent",
                                         createModifiedCallback(self))
        [self._arrayselection.AddArray(k) for k in TAG_RANGES.keys()]
        self._arrayselection.EnableAllArrays()

        self.Modified()


    #### Following methods are for composite datasets ####

    def use_composite(self):
        """True if input dataset is a composite dataset"""
        return self.__composite

    def set_block_filenames(self, n):
        """Gets a list of filenames based on user input filename and creates a
        numbered list of filenames for the reader to save out. Assumes the
        filename has an extension set already.
        """
        number = n
        count = 0
        while number > 0:
            number = number // 10
            count = count + 1
        count = '%d' % count
        identifier = '_%.' + count + 'd'
        blocknum = [identifier % i for i in range(n)]
        # Check the file extension:
        ext = self.get_file_name().split('.')[-1]
        basename = self.get_file_name().replace('.%s' % ext, '')
        self.__blockfilenames = [
            basename + '%s.%s' % (blocknum[i], ext) for i in range(n)
        ]
        return self.__blockfilenames

    def get_block_filename(self, idx):
        """Get the filename for a specific block if composite dataset."""
        return self.__blockfilenames[idx]

    def RequestData(self, request, inInfo, outInfo):
        """Subclasses must implement a ``PerformWriteOut`` method that takes an
        input data object and a filename. This method will automatically handle
        composite data sets.
        """
        inp = self.GetInputData(inInfo, 0, 0)
        if isinstance(inp, vtk.vtkMultiBlockDataSet):
            print("instance is a vtkmultiblockdataset")
            self.__composite = True
        # Handle composite datasets. NOTE: This only handles vtkMultiBlockDataSet
        if self.__composite:
            num = inp.GetNumberOfBlocks()
            print("number of blocks = ",num)
            self.set_block_filenames(num)
            for i in range(num):
                data = inp.GetBlock(i)
                #print("data = ",data)
                name = inp.GetMetaData(i).Get(vtk.vtkCompositeDataSet.NAME())
                print("metadata block name = ",name)
                print("input type = ",self.InputType)
                if data.IsTypeOf(self.InputType):
                    self.PerformWriteOut(data, self.get_block_filename(i), name)
                else:
                    print("input block not of correct type")
                    #warnings.warn(
                    #    'Input block %d of type(%s) not saveable by writer.'
                    #    % (i, type(data))
                    #)
        # Handle single input dataset
        else:
            self.PerformWriteOut(inp, self.get_file_name(), None)
        return 1


###############################################################################
## Now lets use ``WriterBase`` to make a writer algorithm that ParaView can use


class WriteCellCenterData(WriterBase):
    """This writer will save a file of the XYZ points for an input dataset's
    cell centers and its cell data. Use in tandom with ParaView's native CSV
    writer which saves the PointData. This class was originally
    implemented in `PVGeo`_ by `Bane Sullivan`_.

    .. _PVGeo: http://pvgeo.org
    .. _Bane Sullivan: http://banesullivan.com
    """

    def __init__(self):
        # nijso changed from vtkdataset to vtkmultiblockdataset
        WriterBase.__init__(self, inputType='vtkMultiBlockDataSet')

    def PerformWriteOut(self, inp, filename, object_name):
        print("performwriteout, block =",object_name)

        # number of blocks should be 2 (internal and boundary)
        num = inp.GetNumberOfBlocks()
        print("number of blocks = ",num)

        cel = inp.GetNumberOfCells()
        print("number of cells = ",cel)

        pnt = inp.GetNumberOfPoints()
        print("number of points = ",pnt)

        # first, get the internal part
        dataInternal = inp.GetBlock(0)
        #print("name internal = ",inp.GetName(0))
        name = inp.GetMetaData(0).Get(vtk.vtkCompositeDataSet.NAME())
        print("metadata block name = ",name)
        dataBoundary = inp.GetBlock(1)
        #print("name Boundary = ",inp.GetName(1))
        name = inp.GetMetaData(1).Get(vtk.vtkCompositeDataSet.NAME())
        print("metadata block name = ",name)
        numMarkers = dataBoundary.GetNumberOfBlocks()
        print("number of boundary blocks = ",numMarkers)
        markernames = []
        for marker in range(numMarkers):
          name = dataBoundary.GetMetaData(marker).Get(vtk.vtkCompositeDataSet.NAME())
          print("Boundary block name = ",name)
          markernames.append(name)


        # ### ### #

        markerconnectivities=[]
        print("marker connectivities = ",markerconnectivities)

        # lets try something. First we add the connectivity to one of the boundary markers
        for imarker in range(numMarkers):
          dataInlet = dataBoundary.GetBlock(imarker)
          #print("data marker = ",dataInlet)
          #print("datainlet = ",dir(dataInlet))
          #dataInlet.Update()
          #print("update =",dataInlet)
          connect = vtk.vtkConnectivityFilter()
          #input needs to be a vtkAlgorithmoutput
          #connect.SetInputConnection(dataInlet)
          ###   ###
          connect.SetInputData(dataInlet)
          ###   ###
          connect.SetExtractionModeToAllRegions()
          connect.ColorRegionsOn()
          connect.Update()
          #print("connect = ",connect)

          newdataInlet = vtkUnstructuredGrid()
          newDataInlet = connect.GetUnstructuredGridOutput()
          #print("newinlet = ",newDataInlet)
          #print("newinlet = ",dir(newDataInlet))


          # add a calculator to add the RegionId as a field
          calc = vtk.vtkArrayCalculator()
          #calc.SetInput(newDataInlet)
          calc.SetInputData(newDataInlet)
          calc.SetAttributeTypeToPointData()
          # do not ignore datasets where field data array is not present
          calc.IgnoreMissingArraysOff()
          calc.AddScalarArrayName("RegionId")
          calc.SetFunction("RegionId+1")
          #calc.AddScalarArrayName("Pressure")
          #calc.SetFunction("Pressure+1.0")

          calc.SetResultArrayName("newRegionID")

          print("### array names = ",calc.GetScalarArrayNames())
          #print("### array pressure valid = ",calc.CheckValidVariableName("Pressure"))
          print("### array nr = ",calc.GetNumberOfScalarArrays())
          print("### variable names = ",calc.GetScalarVariableNames())

          calc.ReplaceInvalidValuesOn()
          calc.SetReplacementValue(0.0)
          calc.Update() 
          print(calc.GetOutput())

          print("newcalc = ",calc)
          #print("newcalc = ",dir(calc))

          newdataInlet2 = vtkUnstructuredGrid()
          newDataInlet2 = calc.GetUnstructuredGridOutput()
          print("newinlet2 = ",newDataInlet2)
          #print("newinlet2 = ",dir(newDataInlet2))


          # ### we now have the regionID added to the field arrays. Now we need ResampleWithDataset ### #
          # databoundary : newDataInlet2 (vtkUnstructuredgrid)
          # datainternal : dataInternal (multiblockdataset)
          resampler = vtk.vtkResampleWithDataSet()
          # input 
          # source (samples from source to input)
          # output gets type of input
          resampler.SetSourceData(newDataInlet2)
          resampler.SetInputData(dataInternal)
          #resampler.SetSourceConnection(dataInternal)
          resampler.Update()
          print("resampler = ",resampler)
          #print("resampler = ",dir(resampler))
          mboutput =vtkMultiBlockDataSet()
          mboutput = resampler.GetOutput()
          print("mboutput = ",mboutput)
          #print("mboutput = ",dir(mboutput))


          # ### we now have an array newRegionID, let's check what is inside.

          # 1. get the unstructured grid (assuming we have one single internal block)
          mb1 = mboutput.GetBlock(0)
          #print(mb1)
          #print(dir(mb1))
          #print("mbgrid number of cells = ",mb1.GetNumberOfCells())
          #print("mbgrid number of points = ",mb1.GetNumberOfPoints())
          field1 = mb1.GetPointData()
          #print(field1)
          #print(dir(field1))
          #print(field1.GetNumberOfArrays())
          field2 = field1.GetArray('newRegionID')

          #field1 = mb1.GetData('Pressure')
          #print(dir(field2))
          # the points on the marker have value "1"
          # so we now have the global point IDs on the marker, but we need the cell connectivity (actually the edges)
          #print("number of values = ",field2.GetNumberOfValues())

          # contains the indices of the points on the boundary
          setOfPoints= set()
          setOfCells= set()
          # contains a list of connectivities [edge1, edge2] where the inside of the domain is on the left of the edge
          listOfEdgeConnectivities= []

          # construct the set of points on the marker
          for i in range(mb1.GetNumberOfPoints()):
            print(i," , ",field2.GetValue(i))
            # should add to the set if the value is 1.0
            if field2.GetValue(i)>0.5:
                setOfPoints.add(i)

          print("set of points on the boundary of ",markernames[imarker] ," is: ",setOfPoints)

          # we now construct a list of all cells containing the points on the boundary
          # In that way we can construct the connectivity 
          allCellIdsList = vtk.vtkIdList()

          # loop over the set of points
          for pointId in setOfPoints:
              cellIdsList = vtk.vtkIdList()
              # not sure how efficient this is...
              #ug.GetPointCells(pointId, cellIdsList)
              mb1.GetPointCells(pointId, cellIdsList)

              # loop over all cells that contain the point
              for i in range(cellIdsList.GetNumberOfIds()):
                  allCellIdsList.InsertNextId(cellIdsList.GetId(i))
                  setOfCells.add(cellIdsList.GetId(i))

          print("list of all cells containing the points:",setOfCells)
          #print("list of all cells containing the points:",allCellIdsList)
          print(dir(allCellIdsList))

          # now get the edge connectivity
          cellIdList = vtk.vtkIdList()
          ugrid = dataInternal.GetBlock(0)
          celldata = ugrid.GetCells()
          print("celldata = ",dir(celldata))

          singleMarkerConnectivities=[]
          # ### we loop over all cells ### #
          for cellID in setOfCells:
            print("value = ",cellID)  
            print("value = ",type(cellID))  
            celldata.GetCellAtId(cellID,cellIdList)
            print("number of ids = ",cellIdList.GetNumberOfIds())
            print("celltype = ",ugrid.GetCellType(cellID))

            # note that a single cell can have more than one edge on the marker
            edgeID=[]
            # loop over all points of the cell

            singleCellIdList=[]
            for j in range(cellIdList.GetNumberOfIds()):
              singleCellIdList.append(cellIdList.GetId(j))
            print("single cell id list = ",singleCellIdList)  
          
            #for j in range(cellIdList.GetNumberOfIds()):
            for j in range(0,len(singleCellIdList)):
              print("__j =",j," , id=",singleCellIdList[j])

              # check if point of the cell is in list of points on the edge
              if (singleCellIdList[j] in setOfPoints):
                print("on the boundary nr ",j," , id=",singleCellIdList[j])
                if (singleCellIdList[j-1] in setOfPoints):
                  # make sure that the ordering is correct for the influx/outflux  
                  print("complete edge on the boundary =",j," , id=",singleCellIdList[j]," ",singleCellIdList[j-1])
                  singleMarkerConnectivities.append([singleCellIdList[j],singleCellIdList[j-1]])
                  print("single marker connectivities = ",singleMarkerConnectivities)


          markerconnectivities.append(singleMarkerConnectivities)
          print("marker connectivities = ",markerconnectivities)
          ####### end of boundary marker test #####
        


        #print("### data = ",dataInternal)
        #get the unstructured data part (just one level deeper)
        ugrid = dataInternal.GetBlock(0)
        #print("### data 0= ",ugrid)
        #print("###   dir of a unstructuredgrid:   ###")
        #print(dir(ugrid))

        #name = data.GetMetaData(0).Get(vtk.vtkCompositeDataSet.NAME())
        #print("metadata block name = ",name)
        #print("input type = ",self.InputType)

        cellnr = ugrid.GetNumberOfCells()
        print("ugrid number of cells = ",cellnr)
        pntnr = ugrid.GetNumberOfPoints()
        print("ugrid number of points = ",pntnr)


        # get the cell type (id) [WORKING]
        #for i in range(cellnr):
        #    cel = ugrid.GetCellType(i)
        #    print("celltype = ",cel)
       
        #pnts = vtk.vtkIdList()
        # get the cell points (id)
        #for i in range(cellnr):
            # we have to get the number of points (connectivity entries)
            #data0.GetCellPoints(i,pnts)
            #print("cellpoints = ",pnts)
            #print(pnts(0)," ",pnts(1))

        # get the cell type (id)
        #for i in range(cellnr):
        # this is a vtkcellarray:
        celldata = ugrid.GetCells()
        print("cell = ",celldata)
        print("###   dir of a vtkcellarray:   ###")
        #print(dir(celldata))
        print("cellsize = ",celldata.GetNumberOfCells())
        print(celldata.GetConnectivityArray())
        #for c in range(cel.GetNumberOfCells()):
        #    print(cel.GetCell(c))







        # ### this gets the connectivity ### # 
        print("###   get the connectivity   ###")
        cellIdList = vtk.vtkIdList()

        for i in range(celldata.GetNumberOfCells()):
          celldata.GetCellAtId(i,cellIdList)
          print("number of ids = ",cellIdList.GetNumberOfIds())
          print("celltype = ",)
          for j in range(cellIdList.GetNumberOfIds()):
            print(cellIdList.GetId(j))


        # ### this gets the points ### # 
        print("###   get the points   ###")
        pntdata = ugrid.GetPoints()
        #for p in range(pntdata.GetNumberOfPoints()):
        #    print(pntdata.GetPoint(p))
       



       # ######### #################### ########## #
       # ######### WRITE THE FILE       ########## #
       # ######### #################### ########## #


        with open(filename, 'w') as f:
          line="NDIME= 2\n"  
          f.write(line)
          line="NELEM= " + str(cellnr) + "\n"
          f.write(line)

          # add the connectivities
          pts = vtk.vtkIdList()
          for i in range(celldata.GetNumberOfCells()):
            celldata.GetCellAtId(i,pts)

            #print("number of ids = ",pts.GetNumberOfIds())
            #print("celltype = ",)
            nrnodes = pts.GetNumberOfIds()
            # triangle
            if (nrnodes==3):
              line="5 " + str(pts.GetId(0)) + " " + str(pts.GetId(1)) + " " + str(pts.GetId(2)) + " " + str(i) + "\n"
            # quadrilateral
            if (nrnodes==4):
              line="9 " + str(pts.GetId(0)) + " " + str(pts.GetId(1)) + " " + str(pts.GetId(2)) + " " + str(pts.GetId(3)) + " " + str(i) + "\n"
            f.write(line)    


          ###   add the points   ###
          line="NPOIN= " + str(pntnr) + "\n"
          f.write(line)
         
          # global point coordinates
          for p in range(pntdata.GetNumberOfPoints()):
            pointcoord= pntdata.GetPoint(p)
            print("pointcoord=",pointcoord)
            #print(dir(pointcoord))
            # ##### 3D #####
            #line = str(pointcoord[0]) + " " + str(pointcoord[1]) + " " + str(pointcoord[2]) + " " + str(p) + "\n"
            # ##### 2D #####
            line = str(pointcoord[0]) + " " + str(pointcoord[1]) + " " + str(p) + "\n"
            f.write(line)    


          # ############################################## #
          # add the markers
          # ############################################## #
          numMarkers = dataBoundary.GetNumberOfBlocks()
          print("number of boundary blocks = ",numMarkers)
          line="NMARK= " + str(numMarkers) + "\n"
          f.write(line)

          # ############################################## #
          # loop over the boundary markers 
          # ############################################## #
          for iMarker in range(numMarkers):
            name = dataBoundary.GetMetaData(iMarker).Get(vtk.vtkCompositeDataSet.NAME())
            # 1. write the name of the marker
            line="MARKER_TAG= " + str(name) + "\n"
            f.write(line)
            # 2. write the number of elements of the marker
            markerData = dataBoundary.GetBlock(iMarker)
            marker_elems = markerData.GetNumberOfCells()
            line="MARKER_ELEMS= " + str(marker_elems) + "\n"
            f.write(line)
            # 3. write the type and the connectivity 
            print("dir boundary marker = ",dir(markerData))
            celldata = markerData.GetCells()

            # loop over number of cells on the marker 
            for i in range(celldata.GetNumberOfCells()):
              # write the marker type: "3" = 2D line segment 
              line="3 " + str(markerconnectivities[iMarker][i][0]) + " " + str(markerconnectivities[iMarker][i][1]) + "\n"
              f.write(line)
 

            pts = vtk.vtkIdList()
            # loop over number of cells 
            for i in range(celldata.GetNumberOfCells()):
              celldata.GetCellAtId(i,pts)
              print("number of ids = ",pts.GetNumberOfIds())
              print("celltype = ",)
              for j in range(pts.GetNumberOfIds()):
                print(pts.GetId(j))
                #print(pts)
                #print(dir(pts))

          

          #MARKER_ELEMS= 


        #if data.IsA("vtkCompositeDataSet"):
        #    print("we have a vtkCompositedataset")
        #    ustruct = vtk.vtkUnstructuredGrid()
        #    #from paraview.modules.vtkPVVTKExtensionsMisc import vtkMergeBlocks
        #    #from vtk.numpy_interface import dataset_adapter as dsa
        #    #mergeFilter = vtkMergeBlocks()
        #    #mergeFilter.SetInputData(data.VTKObject)
        #    #mergeFilter.Update()
        #    #ustruct = dsa.WrapDataObject(mergeFilter.GetOutput())
        #    #ustruct = MergeBlocks(data)
        #    ustruct = vtkUnstructuredGrid.SafeDownCast(inp.GetBlock(0))
        #    print(ustruct)
        #    print("mergeblock=",dir(ustruct))


        # cell connectivity



        # then, get the boundary part (markers)
        #data = inp.GetBlock(1)
        #print("### data boundary = ",data)
        #name = data.GetMetaData(0).Get(vtk.vtkCompositeDataSet.NAME())
        #print("metadata block name = ",name)
        #print("input type = ",self.InputType)
        #num = data.GetNumberOfBlocks()
        #print("number of blocks = ",num)
        #cel = data.GetNumberOfCells()
        #print("number of cells = ",cel)
        #pnt = data.GetNumberOfPoints()
        #print("number of points = ",pnt)
        #print(dir(data))
        #print("metadata=",data.GetMetaData)
   


        ## Find cell centers
        #filt = vtk.vtkCellCenters()
        #filt.SetInputDataObject(inp)
        #filt.Update()

        #print("filt = ",dir(filt))

        #centers = dsa.WrapDataObject(filt.GetOutput(0)).Points
        ## Get CellData
        #wpdi = dsa.WrapDataObject(inp)
        #celldata = wpdi.CellData
        #keys = celldata.keys()
        ## Save out using numpy
        ##arr = np.zeros((len(centers), 3 + len(keys)))
        #arr[:, 0:3] = centers
        #for i, name in enumerate(keys):
        #    arr[:, i + 3] = celldata[name]
        ## Now write out the data
        ## Clean data titles to make sure they do not contain the delimiter
        #repl = '_' if self.__delimiter != '_' else '-'
        #for i, name in enumerate(keys):
        #    keys[i] = name.replace(self.__delimiter, repl)
        #header = ('%s' % self.__delimiter).join(['X', 'Y', 'Z'] + keys)
        #np.savetxt(
        #    filename,
        #    arr,
        #    header=header,
        #    delimiter=self.__delimiter,
        #    fmt=self.get_format(),
        #    comments='',
        #)
        ## Success for pipeline
        return 1

# ############################################################### #
@smproxy.writer( extensions="su2", file_description="SU2 mesh writing", support_reload=False)
@smproperty.input(name="Input", port_index=0)
@smdomain.datatype(dataTypes=["vtkCompositeDataSet"], composite_data_supported=True)
class PVWriteCellCenterData(WriteCellCenterData):
    """The ``WriteCellCenterData`` class wrapped for use as a plugin in ParaView.
    Be sure that the ``composite_data_supported`` flag is set to ``True``.
    """

    def __init__(self):
        WriteCellCenterData.__init__(self)
        self._arrayselection.AddObserver("ModifiedEvent",
                                         createModifiedCallback(self))
        self._arrayselection.AddArray("one")
        self._arrayselection.AddArray("two")
        self._arrayselection.AddArray("three")
        self._arrayselection.EnableAllArrays()

        WriteCellCenterData._arrayselection.AddObserver("ModifiedEvent",
                                         createModifiedCallback(self))
        WriteCellCenterData._arrayselection.AddArray("one")
        WriteCellCenterData._arrayselection.AddArray("two")
        WriteCellCenterData._arrayselection.AddArray("three")
        WriteCellCenterData._arrayselection.EnableAllArrays()

    @smproperty.stringvector(name="FileName", panel_visibility="never")
    @smdomain.filelist()
    def SetFileName(self, filename):
        """Specify filename for the su2 mesh to write."""
        WriteCellCenterData.SetFileName(self, filename)

    @smproperty.stringvector(name="Format", default_values='%.9e')
    def set_format(self, fmt):
        """Use to set the ASCII format for the writer default is ``'%.9e'``"""
        WriteCellCenterData.set_format(self, fmt)

    @smproperty.stringvector(name="Delimiter", default_values=',')
    def set_delimiter(self, deli):
        """The string delimiter to use"""
        WriteCellCenterData.set_delimiter(self, deli)


    @smproperty.stringvector(name="string", default_values='bla')
    def set_string(self, stringval):
        """The string value to use"""
        WriteCellCenterData.set_string(self, stringval)

   
    @smproperty.dataarrayselection(name="Tag types")
    def GetDataArraySelection(self):
        return self._get_array_selection()
