from __future__ import print_function, absolute_import, division
from . import _MiddleMan
import f90wrap.runtime
import logging
import numpy
import warnings
import weakref

class Middleman(f90wrap.runtime.FortranModule):
    """
    Module middleman
    Defined at MiddleMan.f95 lines 1-1147
    """
    @f90wrap.runtime.register_class("MiddleMan.FortranData")
    class FortranData(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=fortrandata)
        Defined at MiddleMan.f95 lines 4-37
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for fortrandata
            
            self = Fortrandata()
            Defined at MiddleMan.f95 lines 4-37
            
            Returns
            -------
            this : Fortrandata
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = _MiddleMan.f90wrap_middleman__fortrandata_initialise()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_MiddleMan, "f90wrap_middleman__fortrandata_finalise")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def startrigidity(self):
            """
            Element startrigidity ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 6
            """
            return _MiddleMan.f90wrap_fortrandata__get__startrigidity(self._handle)
        
        @startrigidity.setter
        def startrigidity(self, startrigidity):
            _MiddleMan.f90wrap_fortrandata__set__startrigidity(self._handle, startrigidity)
        
        @property
        def endrigidity(self):
            """
            Element endrigidity ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 7
            """
            return _MiddleMan.f90wrap_fortrandata__get__endrigidity(self._handle)
        
        @endrigidity.setter
        def endrigidity(self, endrigidity):
            _MiddleMan.f90wrap_fortrandata__set__endrigidity(self._handle, endrigidity)
        
        @property
        def rigiditystep(self):
            """
            Element rigiditystep ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 8
            """
            return _MiddleMan.f90wrap_fortrandata__get__rigiditystep(self._handle)
        
        @rigiditystep.setter
        def rigiditystep(self, rigiditystep):
            _MiddleMan.f90wrap_fortrandata__set__rigiditystep(self._handle, rigiditystep)
        
        @property
        def gyropercent(self):
            """
            Element gyropercent ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 9
            """
            return _MiddleMan.f90wrap_fortrandata__get__gyropercent(self._handle)
        
        @gyropercent.setter
        def gyropercent(self, gyropercent):
            _MiddleMan.f90wrap_fortrandata__set__gyropercent(self._handle, gyropercent)
        
        @property
        def fixedstepsize(self):
            """
            Element fixedstepsize ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 10
            """
            return _MiddleMan.f90wrap_fortrandata__get__fixedstepsize(self._handle)
        
        @fixedstepsize.setter
        def fixedstepsize(self, fixedstepsize):
            _MiddleMan.f90wrap_fortrandata__set__fixedstepsize(self._handle, fixedstepsize)
        
        @property
        def sphere(self):
            """
            Element sphere ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 11
            """
            return _MiddleMan.f90wrap_fortrandata__get__sphere(self._handle)
        
        @sphere.setter
        def sphere(self, sphere):
            _MiddleMan.f90wrap_fortrandata__set__sphere(self._handle, sphere)
        
        @property
        def trapdist(self):
            """
            Element trapdist ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 12
            """
            return _MiddleMan.f90wrap_fortrandata__get__trapdist(self._handle)
        
        @trapdist.setter
        def trapdist(self, trapdist):
            _MiddleMan.f90wrap_fortrandata__set__trapdist(self._handle, trapdist)
        
        @property
        def berr(self):
            """
            Element berr ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 13
            """
            return _MiddleMan.f90wrap_fortrandata__get__berr(self._handle)
        
        @berr.setter
        def berr(self, berr):
            _MiddleMan.f90wrap_fortrandata__set__berr(self._handle, berr)
        
        @property
        def transmissionrres(self):
            """
            Element transmissionrres ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 14
            """
            return _MiddleMan.f90wrap_fortrandata__get__transmissionrres(self._handle)
        
        @transmissionrres.setter
        def transmissionrres(self, transmissionrres):
            _MiddleMan.f90wrap_fortrandata__set__transmissionrres(self._handle, \
                transmissionrres)
        
        @property
        def gaussianlength(self):
            """
            Element gaussianlength ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 15
            """
            return _MiddleMan.f90wrap_fortrandata__get__gaussianlength(self._handle)
        
        @gaussianlength.setter
        def gaussianlength(self, gaussianlength):
            _MiddleMan.f90wrap_fortrandata__set__gaussianlength(self._handle, \
                gaussianlength)
        
        @property
        def positionin(self):
            """
            Element positionin ftype=real  pytype=float array
            Defined at MiddleMan.f95 line 16
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_fortrandata__array__positionin(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            positionin = self._arrays.get(array_hash)
            if positionin is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if positionin.ctypes.data != array_handle:
                    positionin = None
            if positionin is None:
                try:
                    positionin = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_fortrandata__array__positionin)
                except TypeError:
                    positionin = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = positionin
            return positionin
        
        @positionin.setter
        def positionin(self, positionin):
            self.positionin[...] = positionin
        
        @property
        def date(self):
            """
            Element date ftype=real  pytype=float array
            Defined at MiddleMan.f95 line 17
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_fortrandata__array__date(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            date = self._arrays.get(array_hash)
            if date is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if date.ctypes.data != array_handle:
                    date = None
            if date is None:
                try:
                    date = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_fortrandata__array__date)
                except TypeError:
                    date = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = date
            return date
        
        @date.setter
        def date(self, date):
            self.date[...] = date
        
        @property
        def wind(self):
            """
            Element wind ftype=real  pytype=float array
            Defined at MiddleMan.f95 line 18
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_fortrandata__array__wind(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            wind = self._arrays.get(array_hash)
            if wind is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if wind.ctypes.data != array_handle:
                    wind = None
            if wind is None:
                try:
                    wind = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_fortrandata__array__wind)
                except TypeError:
                    wind = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = wind
            return wind
        
        @wind.setter
        def wind(self, wind):
            self.wind[...] = wind
        
        @property
        def end(self):
            """
            Element end ftype=real  pytype=float array
            Defined at MiddleMan.f95 line 19
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_fortrandata__array__end(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            end = self._arrays.get(array_hash)
            if end is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if end.ctypes.data != array_handle:
                    end = None
            if end is None:
                try:
                    end = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_fortrandata__array__end)
                except TypeError:
                    end = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = end
            return end
        
        @end.setter
        def end(self, end):
            self.end[...] = end
        
        @property
        def intmode(self):
            """
            Element intmode ftype=integer(8) pytype=int32
            Defined at MiddleMan.f95 line 20
            """
            return _MiddleMan.f90wrap_fortrandata__get__intmode(self._handle)
        
        @intmode.setter
        def intmode(self, intmode):
            _MiddleMan.f90wrap_fortrandata__set__intmode(self._handle, intmode)
        
        @property
        def atomicnumber(self):
            """
            Element atomicnumber ftype=integer(8) pytype=int32
            Defined at MiddleMan.f95 line 21
            """
            return _MiddleMan.f90wrap_fortrandata__get__atomicnumber(self._handle)
        
        @atomicnumber.setter
        def atomicnumber(self, atomicnumber):
            _MiddleMan.f90wrap_fortrandata__set__atomicnumber(self._handle, atomicnumber)
        
        @property
        def anti(self):
            """
            Element anti ftype=integer(8) pytype=int32
            Defined at MiddleMan.f95 line 22
            """
            return _MiddleMan.f90wrap_fortrandata__get__anti(self._handle)
        
        @anti.setter
        def anti(self, anti):
            _MiddleMan.f90wrap_fortrandata__set__anti(self._handle, anti)
        
        @property
        def iopt(self):
            """
            Element iopt ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 23
            """
            return _MiddleMan.f90wrap_fortrandata__get__iopt(self._handle)
        
        @iopt.setter
        def iopt(self, iopt):
            _MiddleMan.f90wrap_fortrandata__set__iopt(self._handle, iopt)
        
        @property
        def mode(self):
            """
            Element mode ftype=integer(4) pytype=int array
            Defined at MiddleMan.f95 line 24
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_fortrandata__array__mode(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            mode = self._arrays.get(array_hash)
            if mode is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if mode.ctypes.data != array_handle:
                    mode = None
            if mode is None:
                try:
                    mode = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_fortrandata__array__mode)
                except TypeError:
                    mode = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = mode
            return mode
        
        @mode.setter
        def mode(self, mode):
            self.mode[...] = mode
        
        @property
        def pause(self):
            """
            Element pause ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 25
            """
            return _MiddleMan.f90wrap_fortrandata__get__pause(self._handle)
        
        @pause.setter
        def pause(self, pause):
            _MiddleMan.f90wrap_fortrandata__set__pause(self._handle, pause)
        
        @property
        def rcomputation(self):
            """
            Element rcomputation ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 26
            """
            return _MiddleMan.f90wrap_fortrandata__get__rcomputation(self._handle)
        
        @rcomputation.setter
        def rcomputation(self, rcomputation):
            _MiddleMan.f90wrap_fortrandata__set__rcomputation(self._handle, rcomputation)
        
        @property
        def scanchoice(self):
            """
            Element scanchoice ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 27
            """
            return _MiddleMan.f90wrap_fortrandata__get__scanchoice(self._handle)
        
        @scanchoice.setter
        def scanchoice(self, scanchoice):
            _MiddleMan.f90wrap_fortrandata__set__scanchoice(self._handle, scanchoice)
        
        @property
        def fortranthreads(self):
            """
            Element fortranthreads ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 28
            """
            return _MiddleMan.f90wrap_fortrandata__get__fortranthreads(self._handle)
        
        @fortranthreads.setter
        def fortranthreads(self, fortranthreads):
            _MiddleMan.f90wrap_fortrandata__set__fortranthreads(self._handle, \
                fortranthreads)
        
        @property
        def n(self):
            """
            Element n ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 29
            """
            return _MiddleMan.f90wrap_fortrandata__get__n(self._handle)
        
        @n.setter
        def n(self, n):
            _MiddleMan.f90wrap_fortrandata__set__n(self._handle, n)
        
        @property
        def transmissionsamples(self):
            """
            Element transmissionsamples ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 30
            """
            return _MiddleMan.f90wrap_fortrandata__get__transmissionsamples(self._handle)
        
        @transmissionsamples.setter
        def transmissionsamples(self, transmissionsamples):
            _MiddleMan.f90wrap_fortrandata__set__transmissionsamples(self._handle, \
                transmissionsamples)
        
        @property
        def maxdegree(self):
            """
            Element maxdegree ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 31
            """
            return _MiddleMan.f90wrap_fortrandata__get__maxdegree(self._handle)
        
        @maxdegree.setter
        def maxdegree(self, maxdegree):
            _MiddleMan.f90wrap_fortrandata__set__maxdegree(self._handle, maxdegree)
        
        @property
        def coordsystem(self):
            """
            Element coordsystem ftype=character(len=3) pytype=str
            Defined at MiddleMan.f95 line 32
            """
            return _MiddleMan.f90wrap_fortrandata__get__coordsystem(self._handle)
        
        @coordsystem.setter
        def coordsystem(self, coordsystem):
            _MiddleMan.f90wrap_fortrandata__set__coordsystem(self._handle, coordsystem)
        
        @property
        def mhdcoordsys(self):
            """
            Element mhdcoordsys ftype=character(len=3) pytype=str
            Defined at MiddleMan.f95 line 33
            """
            return _MiddleMan.f90wrap_fortrandata__get__mhdcoordsys(self._handle)
        
        @mhdcoordsys.setter
        def mhdcoordsys(self, mhdcoordsys):
            _MiddleMan.f90wrap_fortrandata__set__mhdcoordsys(self._handle, mhdcoordsys)
        
        @property
        def inputcoord(self):
            """
            Element inputcoord ftype=character(len=3) pytype=str
            Defined at MiddleMan.f95 line 34
            """
            return _MiddleMan.f90wrap_fortrandata__get__inputcoord(self._handle)
        
        @inputcoord.setter
        def inputcoord(self, inputcoord):
            _MiddleMan.f90wrap_fortrandata__set__inputcoord(self._handle, inputcoord)
        
        @property
        def adapt(self):
            """
            Element adapt ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 35
            """
            return _MiddleMan.f90wrap_fortrandata__get__adapt(self._handle)
        
        @adapt.setter
        def adapt(self, adapt):
            _MiddleMan.f90wrap_fortrandata__set__adapt(self._handle, adapt)
        
        @property
        def totalbetacheck(self):
            """
            Element totalbetacheck ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 36
            """
            return _MiddleMan.f90wrap_fortrandata__get__totalbetacheck(self._handle)
        
        @totalbetacheck.setter
        def totalbetacheck(self, totalbetacheck):
            _MiddleMan.f90wrap_fortrandata__set__totalbetacheck(self._handle, \
                totalbetacheck)
        
        @property
        def trapdistcheck(self):
            """
            Element trapdistcheck ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 37
            """
            return _MiddleMan.f90wrap_fortrandata__get__trapdistcheck(self._handle)
        
        @trapdistcheck.setter
        def trapdistcheck(self, trapdistcheck):
            _MiddleMan.f90wrap_fortrandata__set__trapdistcheck(self._handle, trapdistcheck)
        
        def __str__(self):
            ret = ['<fortrandata>{\n']
            ret.append('    startrigidity : ')
            ret.append(repr(self.startrigidity))
            ret.append(',\n    endrigidity : ')
            ret.append(repr(self.endrigidity))
            ret.append(',\n    rigiditystep : ')
            ret.append(repr(self.rigiditystep))
            ret.append(',\n    gyropercent : ')
            ret.append(repr(self.gyropercent))
            ret.append(',\n    fixedstepsize : ')
            ret.append(repr(self.fixedstepsize))
            ret.append(',\n    sphere : ')
            ret.append(repr(self.sphere))
            ret.append(',\n    trapdist : ')
            ret.append(repr(self.trapdist))
            ret.append(',\n    berr : ')
            ret.append(repr(self.berr))
            ret.append(',\n    transmissionrres : ')
            ret.append(repr(self.transmissionrres))
            ret.append(',\n    gaussianlength : ')
            ret.append(repr(self.gaussianlength))
            ret.append(',\n    positionin : ')
            ret.append(repr(self.positionin))
            ret.append(',\n    date : ')
            ret.append(repr(self.date))
            ret.append(',\n    wind : ')
            ret.append(repr(self.wind))
            ret.append(',\n    end : ')
            ret.append(repr(self.end))
            ret.append(',\n    intmode : ')
            ret.append(repr(self.intmode))
            ret.append(',\n    atomicnumber : ')
            ret.append(repr(self.atomicnumber))
            ret.append(',\n    anti : ')
            ret.append(repr(self.anti))
            ret.append(',\n    iopt : ')
            ret.append(repr(self.iopt))
            ret.append(',\n    mode : ')
            ret.append(repr(self.mode))
            ret.append(',\n    pause : ')
            ret.append(repr(self.pause))
            ret.append(',\n    rcomputation : ')
            ret.append(repr(self.rcomputation))
            ret.append(',\n    scanchoice : ')
            ret.append(repr(self.scanchoice))
            ret.append(',\n    fortranthreads : ')
            ret.append(repr(self.fortranthreads))
            ret.append(',\n    n : ')
            ret.append(repr(self.n))
            ret.append(',\n    transmissionsamples : ')
            ret.append(repr(self.transmissionsamples))
            ret.append(',\n    maxdegree : ')
            ret.append(repr(self.maxdegree))
            ret.append(',\n    coordsystem : ')
            ret.append(repr(self.coordsystem))
            ret.append(',\n    mhdcoordsys : ')
            ret.append(repr(self.mhdcoordsys))
            ret.append(',\n    inputcoord : ')
            ret.append(repr(self.inputcoord))
            ret.append(',\n    adapt : ')
            ret.append(repr(self.adapt))
            ret.append(',\n    totalbetacheck : ')
            ret.append(repr(self.totalbetacheck))
            ret.append(',\n    trapdistcheck : ')
            ret.append(repr(self.trapdistcheck))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @f90wrap.runtime.register_class("MiddleMan.ParticleData")
    class ParticleData(f90wrap.runtime.FortranDerivedType):
        """
        Type(name=particledata)
        Defined at MiddleMan.f95 lines 39-94
        """
        def __init__(self, handle=None):
            """
            Automatically generated constructor for particledata
            
            self = Particledata()
            Defined at MiddleMan.f95 lines 39-94
            
            Returns
            -------
            this : Particledata
                Object to be constructed
            
            """
            f90wrap.runtime.FortranDerivedType.__init__(self)
            if isinstance(handle, numpy.ndarray) and handle.ndim == 1 and handle.dtype.num \
                == 5:
                self._handle = handle
                self._alloc = True
            else:
                result = _MiddleMan.f90wrap_middleman__particledata_initialise()
                self._handle = result[0] if isinstance(result, tuple) else result
                self._alloc = True
            self._setup_finalizer()
        
        def _setup_finalizer(self):
            """Set up weak reference destructor to prevent Fortran memory leaks."""
            if self._alloc:
                destructor = getattr(_MiddleMan, "f90wrap_middleman__particledata_finalise")
                self._finalizer = weakref.finalize(self, destructor, self._handle)
        
        @property
        def velocity(self):
            """
            Element velocity ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 41
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__velocity(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            velocity = self._arrays.get(array_hash)
            if velocity is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if velocity.ctypes.data != array_handle:
                    velocity = None
            if velocity is None:
                try:
                    velocity = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__velocity)
                except TypeError:
                    velocity = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = velocity
            return velocity
        
        @velocity.setter
        def velocity(self, velocity):
            self.velocity[...] = velocity
        
        @property
        def geovelocity(self):
            """
            Element geovelocity ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 42
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__geovelocity(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            geovelocity = self._arrays.get(array_hash)
            if geovelocity is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if geovelocity.ctypes.data != array_handle:
                    geovelocity = None
            if geovelocity is None:
                try:
                    geovelocity = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__geovelocity)
                except TypeError:
                    geovelocity = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = geovelocity
            return geovelocity
        
        @geovelocity.setter
        def geovelocity(self, geovelocity):
            self.geovelocity[...] = geovelocity
        
        @property
        def positionarray(self):
            """
            Element positionarray ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 43
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__positionarray(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            positionarray = self._arrays.get(array_hash)
            if positionarray is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if positionarray.ctypes.data != array_handle:
                    positionarray = None
            if positionarray is None:
                try:
                    positionarray = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__positionarray)
                except TypeError:
                    positionarray = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = positionarray
            return positionarray
        
        @positionarray.setter
        def positionarray(self, positionarray):
            self.positionarray[...] = positionarray
        
        @property
        def oldpositionarray(self):
            """
            Element oldpositionarray ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 47
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__oldpositionarray(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            oldpositionarray = self._arrays.get(array_hash)
            if oldpositionarray is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if oldpositionarray.ctypes.data != array_handle:
                    oldpositionarray = None
            if oldpositionarray is None:
                try:
                    oldpositionarray = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__oldpositionarray)
                except TypeError:
                    oldpositionarray = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = oldpositionarray
            return oldpositionarray
        
        @oldpositionarray.setter
        def oldpositionarray(self, oldpositionarray):
            self.oldpositionarray[...] = oldpositionarray
        
        @property
        def velocityarray(self):
            """
            Element velocityarray ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 51
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__velocityarray(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            velocityarray = self._arrays.get(array_hash)
            if velocityarray is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if velocityarray.ctypes.data != array_handle:
                    velocityarray = None
            if velocityarray is None:
                try:
                    velocityarray = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__velocityarray)
                except TypeError:
                    velocityarray = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = velocityarray
            return velocityarray
        
        @velocityarray.setter
        def velocityarray(self, velocityarray):
            self.velocityarray[...] = velocityarray
        
        @property
        def oldvelocityarray(self):
            """
            Element oldvelocityarray ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 54
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__oldvelocityarray(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            oldvelocityarray = self._arrays.get(array_hash)
            if oldvelocityarray is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if oldvelocityarray.ctypes.data != array_handle:
                    oldvelocityarray = None
            if oldvelocityarray is None:
                try:
                    oldvelocityarray = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__oldvelocityarray)
                except TypeError:
                    oldvelocityarray = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = oldvelocityarray
            return oldvelocityarray
        
        @oldvelocityarray.setter
        def oldvelocityarray(self, oldvelocityarray):
            self.oldvelocityarray[...] = oldvelocityarray
        
        @property
        def m(self):
            """
            Element m ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 58
            """
            return _MiddleMan.f90wrap_particledata__get__m(self._handle)
        
        @m.setter
        def m(self, m):
            _MiddleMan.f90wrap_particledata__set__m(self._handle, m)
        
        @property
        def q(self):
            """
            Element q ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 59
            """
            return _MiddleMan.f90wrap_particledata__get__q(self._handle)
        
        @q.setter
        def q(self, q):
            _MiddleMan.f90wrap_particledata__set__q(self._handle, q)
        
        @property
        def z(self):
            """
            Element z ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 60
            """
            return _MiddleMan.f90wrap_particledata__get__z(self._handle)
        
        @z.setter
        def z(self, z):
            _MiddleMan.f90wrap_particledata__set__z(self._handle, z)
        
        @property
        def a(self):
            """
            Element a ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 61
            """
            return _MiddleMan.f90wrap_particledata__get__a(self._handle)
        
        @a.setter
        def a(self, a):
            _MiddleMan.f90wrap_particledata__set__a(self._handle, a)
        
        @property
        def lat(self):
            """
            Element lat ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 62
            """
            return _MiddleMan.f90wrap_particledata__get__lat(self._handle)
        
        @lat.setter
        def lat(self, lat):
            _MiddleMan.f90wrap_particledata__set__lat(self._handle, lat)
        
        @property
        def long_bn(self):
            """
            Element long_bn ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 63
            """
            return _MiddleMan.f90wrap_particledata__get__long_bn(self._handle)
        
        @long_bn.setter
        def long_bn(self, long_bn):
            _MiddleMan.f90wrap_particledata__set__long_bn(self._handle, long_bn)
        
        @property
        def e_0(self):
            """
            Element e_0 ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 64
            """
            return _MiddleMan.f90wrap_particledata__get__e_0(self._handle)
        
        @e_0.setter
        def e_0(self, e_0):
            _MiddleMan.f90wrap_particledata__set__e_0(self._handle, e_0)
        
        @property
        def r(self):
            """
            Element r ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 65
            """
            return _MiddleMan.f90wrap_particledata__get__r(self._handle)
        
        @r.setter
        def r(self, r):
            _MiddleMan.f90wrap_particledata__set__r(self._handle, r)
        
        @property
        def lambda_(self):
            """
            Element lambda_ ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 66
            """
            return _MiddleMan.f90wrap_particledata__get__lambda_(self._handle)
        
        @lambda_.setter
        def lambda_(self, lambda_):
            _MiddleMan.f90wrap_particledata__set__lambda_(self._handle, lambda_)
        
        @property
        def secondtotal(self):
            """
            Element secondtotal ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 68
            """
            return _MiddleMan.f90wrap_particledata__get__secondtotal(self._handle)
        
        @secondtotal.setter
        def secondtotal(self, secondtotal):
            _MiddleMan.f90wrap_particledata__set__secondtotal(self._handle, secondtotal)
        
        @property
        def oldsecondtotal(self):
            """
            Element oldsecondtotal ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 69
            """
            return _MiddleMan.f90wrap_particledata__get__oldsecondtotal(self._handle)
        
        @oldsecondtotal.setter
        def oldsecondtotal(self, oldsecondtotal):
            _MiddleMan.f90wrap_particledata__set__oldsecondtotal(self._handle, \
                oldsecondtotal)
        
        @property
        def timeelapsed(self):
            """
            Element timeelapsed ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 70
            """
            return _MiddleMan.f90wrap_particledata__get__timeelapsed(self._handle)
        
        @timeelapsed.setter
        def timeelapsed(self, timeelapsed):
            _MiddleMan.f90wrap_particledata__set__timeelapsed(self._handle, timeelapsed)
        
        @property
        def h(self):
            """
            Element h ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 72
            """
            return _MiddleMan.f90wrap_particledata__get__h(self._handle)
        
        @h.setter
        def h(self, h):
            _MiddleMan.f90wrap_particledata__set__h(self._handle, h)
        
        @property
        def hold(self):
            """
            Element hold ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 73
            """
            return _MiddleMan.f90wrap_particledata__get__hold(self._handle)
        
        @hold.setter
        def hold(self, hold):
            _MiddleMan.f90wrap_particledata__set__hold(self._handle, hold)
        
        @property
        def lasth(self):
            """
            Element lasth ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 74
            """
            return _MiddleMan.f90wrap_particledata__get__lasth(self._handle)
        
        @lasth.setter
        def lasth(self, lasth):
            _MiddleMan.f90wrap_particledata__set__lasth(self._handle, lasth)
        
        @property
        def firsth(self):
            """
            Element firsth ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 75
            """
            return _MiddleMan.f90wrap_particledata__get__firsth(self._handle)
        
        @firsth.setter
        def firsth(self, firsth):
            _MiddleMan.f90wrap_particledata__set__firsth(self._handle, firsth)
        
        @property
        def maxgyropercent(self):
            """
            Element maxgyropercent ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 76
            """
            return _MiddleMan.f90wrap_particledata__get__maxgyropercent(self._handle)
        
        @maxgyropercent.setter
        def maxgyropercent(self, maxgyropercent):
            _MiddleMan.f90wrap_particledata__set__maxgyropercent(self._handle, \
                maxgyropercent)
        
        @property
        def cachedbfield(self):
            """
            Element cachedbfield ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 79
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__cachedbfield(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            cachedbfield = self._arrays.get(array_hash)
            if cachedbfield is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if cachedbfield.ctypes.data != array_handle:
                    cachedbfield = None
            if cachedbfield is None:
                try:
                    cachedbfield = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__cachedbfield)
                except TypeError:
                    cachedbfield = f90wrap.runtime.direct_c_array(array_type, array_shape, \
                        array_handle)
                self._arrays[array_hash] = cachedbfield
            return cachedbfield
        
        @cachedbfield.setter
        def cachedbfield(self, cachedbfield):
            self.cachedbfield[...] = cachedbfield
        
        @property
        def cachedbfieldvalid(self):
            """
            Element cachedbfieldvalid ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 80
            """
            return _MiddleMan.f90wrap_particledata__get__cachedbfieldvalid(self._handle)
        
        @cachedbfieldvalid.setter
        def cachedbfieldvalid(self, cachedbfieldvalid):
            _MiddleMan.f90wrap_particledata__set__cachedbfieldvalid(self._handle, \
                cachedbfieldvalid)
        
        @property
        def distancetraveled(self):
            """
            Element distancetraveled ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 81
            """
            return _MiddleMan.f90wrap_particledata__get__distancetraveled(self._handle)
        
        @distancetraveled.setter
        def distancetraveled(self, distancetraveled):
            _MiddleMan.f90wrap_particledata__set__distancetraveled(self._handle, \
                distancetraveled)
        
        @property
        def mdp(self):
            """
            Element mdp ftype=real(8) pytype=float array
            Defined at MiddleMan.f95 line 83
            """
            array_ndim, array_type, array_shape, array_handle = \
                _MiddleMan.f90wrap_particledata__array__mdp(self._handle)
            array_hash = hash((array_ndim, array_type, tuple(array_shape), array_handle))
            mdp = self._arrays.get(array_hash)
            if mdp is not None:
                # Validate cached array: check data pointer matches current handle (issue #222)
                # Arrays can be deallocated and reallocated at same address, invalidating cache
                if mdp.ctypes.data != array_handle:
                    mdp = None
            if mdp is None:
                try:
                    mdp = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                            self._handle,
                                            _MiddleMan.f90wrap_particledata__array__mdp)
                except TypeError:
                    mdp = f90wrap.runtime.direct_c_array(array_type, array_shape, array_handle)
                self._arrays[array_hash] = mdp
            return mdp
        
        @mdp.setter
        def mdp(self, mdp):
            self.mdp[...] = mdp
        
        @property
        def betaerror(self):
            """
            Element betaerror ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 84
            """
            return _MiddleMan.f90wrap_particledata__get__betaerror(self._handle)
        
        @betaerror.setter
        def betaerror(self, betaerror):
            _MiddleMan.f90wrap_particledata__set__betaerror(self._handle, betaerror)
        
        @property
        def originalbeta(self):
            """
            Element originalbeta ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 85
            """
            return _MiddleMan.f90wrap_particledata__get__originalbeta(self._handle)
        
        @originalbeta.setter
        def originalbeta(self, originalbeta):
            _MiddleMan.f90wrap_particledata__set__originalbeta(self._handle, originalbeta)
        
        @property
        def currentbeta(self):
            """
            Element currentbeta ftype=real(8) pytype=float32
            Defined at MiddleMan.f95 line 86
            """
            return _MiddleMan.f90wrap_particledata__get__currentbeta(self._handle)
        
        @currentbeta.setter
        def currentbeta(self, currentbeta):
            _MiddleMan.f90wrap_particledata__set__currentbeta(self._handle, currentbeta)
        
        @property
        def finalstep(self):
            """
            Element finalstep ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 88
            """
            return _MiddleMan.f90wrap_particledata__get__finalstep(self._handle)
        
        @finalstep.setter
        def finalstep(self, finalstep):
            _MiddleMan.f90wrap_particledata__set__finalstep(self._handle, finalstep)
        
        @property
        def mindistcheck(self):
            """
            Element mindistcheck ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 89
            """
            return _MiddleMan.f90wrap_particledata__get__mindistcheck(self._handle)
        
        @mindistcheck.setter
        def mindistcheck(self, mindistcheck):
            _MiddleMan.f90wrap_particledata__set__mindistcheck(self._handle, mindistcheck)
        
        @property
        def escaped(self):
            """
            Element escaped ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 90
            """
            return _MiddleMan.f90wrap_particledata__get__escaped(self._handle)
        
        @escaped.setter
        def escaped(self, escaped):
            _MiddleMan.f90wrap_particledata__set__escaped(self._handle, escaped)
        
        @property
        def totalbetachecktrigger(self):
            """
            Element totalbetachecktrigger ftype=logical pytype=bool
            Defined at MiddleMan.f95 line 91
            """
            return _MiddleMan.f90wrap_particledata__get__totalbetachecktrigger(self._handle)
        
        @totalbetachecktrigger.setter
        def totalbetachecktrigger(self, totalbetachecktrigger):
            _MiddleMan.f90wrap_particledata__set__totalbetachecktrigger(self._handle, \
                totalbetachecktrigger)
        
        @property
        def steps(self):
            """
            Element steps ftype=integer(8) pytype=int32
            Defined at MiddleMan.f95 line 92
            """
            return _MiddleMan.f90wrap_particledata__get__steps(self._handle)
        
        @steps.setter
        def steps(self, steps):
            _MiddleMan.f90wrap_particledata__set__steps(self._handle, steps)
        
        @property
        def counter(self):
            """
            Element counter ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 93
            """
            return _MiddleMan.f90wrap_particledata__get__counter(self._handle)
        
        @counter.setter
        def counter(self, counter):
            _MiddleMan.f90wrap_particledata__set__counter(self._handle, counter)
        
        @property
        def termtype(self):
            """
            Element termtype ftype=integer(4) pytype=int32
            Defined at MiddleMan.f95 line 94
            """
            return _MiddleMan.f90wrap_particledata__get__termtype(self._handle)
        
        @termtype.setter
        def termtype(self, termtype):
            _MiddleMan.f90wrap_particledata__set__termtype(self._handle, termtype)
        
        def __str__(self):
            ret = ['<particledata>{\n']
            ret.append('    velocity : ')
            ret.append(repr(self.velocity))
            ret.append(',\n    geovelocity : ')
            ret.append(repr(self.geovelocity))
            ret.append(',\n    positionarray : ')
            ret.append(repr(self.positionarray))
            ret.append(',\n    oldpositionarray : ')
            ret.append(repr(self.oldpositionarray))
            ret.append(',\n    velocityarray : ')
            ret.append(repr(self.velocityarray))
            ret.append(',\n    oldvelocityarray : ')
            ret.append(repr(self.oldvelocityarray))
            ret.append(',\n    m : ')
            ret.append(repr(self.m))
            ret.append(',\n    q : ')
            ret.append(repr(self.q))
            ret.append(',\n    z : ')
            ret.append(repr(self.z))
            ret.append(',\n    a : ')
            ret.append(repr(self.a))
            ret.append(',\n    lat : ')
            ret.append(repr(self.lat))
            ret.append(',\n    long_bn : ')
            ret.append(repr(self.long_bn))
            ret.append(',\n    e_0 : ')
            ret.append(repr(self.e_0))
            ret.append(',\n    r : ')
            ret.append(repr(self.r))
            ret.append(',\n    lambda_ : ')
            ret.append(repr(self.lambda_))
            ret.append(',\n    secondtotal : ')
            ret.append(repr(self.secondtotal))
            ret.append(',\n    oldsecondtotal : ')
            ret.append(repr(self.oldsecondtotal))
            ret.append(',\n    timeelapsed : ')
            ret.append(repr(self.timeelapsed))
            ret.append(',\n    h : ')
            ret.append(repr(self.h))
            ret.append(',\n    hold : ')
            ret.append(repr(self.hold))
            ret.append(',\n    lasth : ')
            ret.append(repr(self.lasth))
            ret.append(',\n    firsth : ')
            ret.append(repr(self.firsth))
            ret.append(',\n    maxgyropercent : ')
            ret.append(repr(self.maxgyropercent))
            ret.append(',\n    cachedbfield : ')
            ret.append(repr(self.cachedbfield))
            ret.append(',\n    cachedbfieldvalid : ')
            ret.append(repr(self.cachedbfieldvalid))
            ret.append(',\n    distancetraveled : ')
            ret.append(repr(self.distancetraveled))
            ret.append(',\n    mdp : ')
            ret.append(repr(self.mdp))
            ret.append(',\n    betaerror : ')
            ret.append(repr(self.betaerror))
            ret.append(',\n    originalbeta : ')
            ret.append(repr(self.originalbeta))
            ret.append(',\n    currentbeta : ')
            ret.append(repr(self.currentbeta))
            ret.append(',\n    finalstep : ')
            ret.append(repr(self.finalstep))
            ret.append(',\n    mindistcheck : ')
            ret.append(repr(self.mindistcheck))
            ret.append(',\n    escaped : ')
            ret.append(repr(self.escaped))
            ret.append(',\n    totalbetachecktrigger : ')
            ret.append(repr(self.totalbetachecktrigger))
            ret.append(',\n    steps : ')
            ret.append(repr(self.steps))
            ret.append(',\n    counter : ')
            ret.append(repr(self.counter))
            ret.append(',\n    termtype : ')
            ret.append(repr(self.termtype))
            ret.append('}')
            return ''.join(ret)
        
        _dt_array_initialisers = []
        
    
    @staticmethod
    def cutoff(self, g8, h8, rigidities, allowed, interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        cutoff(self, g8, h8, rigidities, allowed)
        Defined at MiddleMan.f95 lines 106-229
        
        Parameters
        ----------
        data : Fortrandata
        g8 : float array
        h8 : float array
        rigidities : float array
        allowed : int array
        """
        _MiddleMan.f90wrap_middleman__cutoff(data=self._handle, g8=g8, h8=h8, \
            rigidities=rigidities, allowed=allowed)
    
    @staticmethod
    def cone(self, g8, h8, rigidities, allowed, asymlat, asymlong, \
        interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        cone(self, g8, h8, rigidities, allowed, asymlat, asymlong)
        Defined at MiddleMan.f95 lines 242-379
        
        Parameters
        ----------
        data : Fortrandata
        g8 : float array
        h8 : float array
        rigidities : float array
        allowed : int array
        asymlat : float array
        asymlong : float array
        """
        _MiddleMan.f90wrap_middleman__cone(data=self._handle, g8=g8, h8=h8, \
            rigidities=rigidities, allowed=allowed, asymlat=asymlat, asymlong=asymlong)
    
    @staticmethod
    def trajectory_full(self, g8, h8, rigidity, trajectoryfile, trajectoryfilelen, \
        interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        filter, alat, along = trajectory_full(self, g8, h8, rigidity, trajectoryfile, \
            trajectoryfilelen)
        Defined at MiddleMan.f95 lines 394-535
        
        Parameters
        ----------
        data : Fortrandata
        g8 : float array
        h8 : float array
        rigidity : float32
        trajectoryfile : str
        trajectoryfilelen : int32
        
        Returns
        -------
        filter : int32
        alat : float32
        along : float32
        """
        filter, alat, along = \
            _MiddleMan.f90wrap_middleman__trajectory_full(data=self._handle, g8=g8, \
            h8=h8, rigidity=rigidity, trajectoryfile=trajectoryfile, \
            trajectoryfilelen=trajectoryfilelen)
        return filter, alat, along
    
    @staticmethod
    def trajectory(self, g8, h8, rigidities, rigiditieslen, allowed, asymlat, \
        asymlong, interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        trajectory(self, g8, h8, rigidities, rigiditieslen, allowed, asymlat, asymlong)
        Defined at MiddleMan.f95 lines 539-675
        
        Parameters
        ----------
        data : Fortrandata
        g8 : float array
        h8 : float array
        rigidities : float array
        rigiditieslen : int32
        allowed : int array
        asymlat : float array
        asymlong : float array
        """
        _MiddleMan.f90wrap_middleman__trajectory(data=self._handle, g8=g8, h8=h8, \
            rigidities=rigidities, rigiditieslen=rigiditieslen, allowed=allowed, \
            asymlat=asymlat, asymlong=asymlong)
    
    @staticmethod
    def transmission(self, g8, h8, rigidities, transmissions, interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        transmission(self, g8, h8, rigidities, transmissions)
        Defined at MiddleMan.f95 lines 686-830
        
        Parameters
        ----------
        data : Fortrandata
        g8 : float array
        h8 : float array
        rigidities : float array
        transmissions : float array
        """
        _MiddleMan.f90wrap_middleman__transmission(data=self._handle, g8=g8, h8=h8, \
            rigidities=rigidities, transmissions=transmissions)
    
    @staticmethod
    def magstrength(pin, data, coordin, coordout, g8, h8, bfield, \
        interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        magstrength(pin, data, coordin, coordout, g8, h8, bfield)
        Defined at MiddleMan.f95 lines 839-895
        
        Parameters
        ----------
        pin : float array
        data : Fortrandata
        coordin : str
        coordout : str
        g8 : float array
        h8 : float array
        bfield : float array
        """
        _MiddleMan.f90wrap_middleman__magstrength(pin=pin, data=data._handle, \
            coordin=coordin, coordout=coordout, g8=g8, h8=h8, bfield=bfield)
    
    @staticmethod
    def coordtrans(pin, data, coordin, coordout, g8, h8, pout, \
        interface_call=False):
        """
        coordtrans(pin, data, coordin, coordout, g8, h8, pout)
        Defined at MiddleMan.f95 lines 904-929
        
        Parameters
        ----------
        pin : float array
        data : Fortrandata
        coordin : str
        coordout : str
        g8 : float array
        h8 : float array
        pout : float array
        """
        _MiddleMan.f90wrap_middleman__coordtrans(pin=pin, data=data._handle, \
            coordin=coordin, coordout=coordout, g8=g8, h8=h8, pout=pout)
    
    @staticmethod
    def fieldtrace(self, filename, filenamelen, g8, h8, interface_call=False):
        """
        ------------------------------------------------------------------
         Start timing
        ------------------------------------------------------------------
        
        fieldtrace(self, filename, filenamelen, g8, h8)
        Defined at MiddleMan.f95 lines 938-1062
        
        Parameters
        ----------
        data : Fortrandata
        filename : str
        filenamelen : int32
        g8 : float array
        h8 : float array
        """
        _MiddleMan.f90wrap_middleman__fieldtrace(data=self._handle, filename=filename, \
            filenamelen=filenamelen, g8=g8, h8=h8)
    
    @staticmethod
    def mhdstartupsorted(xu, yu, zu, mhdposition_in, mhdb_in, nx_split, ny_split, \
        nz_split, mix, max_bn, miy, may, miz, maz, region_order_in, start_x, end_x, \
        start_y, end_y, start_z, end_z, num_regions, xulen, yulen, zulen, \
        interface_call=False):
        """
        mhdstartupsorted(xu, yu, zu, mhdposition_in, mhdb_in, nx_split, ny_split, \
            nz_split, mix, max_bn, miy, may, miz, maz, region_order_in, start_x, end_x, \
            start_y, end_y, start_z, end_z, num_regions, xulen, yulen, zulen)
        Defined at MiddleMan.f95 lines 1068-1115
        
        Parameters
        ----------
        xu : float array
        yu : float array
        zu : float array
        mhdposition_in : float array
        mhdb_in : float array
        nx_split : int32
        ny_split : int32
        nz_split : int32
        mix : float32
        max_bn : float32
        miy : float32
        may : float32
        miz : float32
        maz : float32
        region_order_in : int array
        start_x : int array
        end_x : int array
        start_y : int array
        end_y : int array
        start_z : int array
        end_z : int array
        num_regions : int32
        xulen : int32
        yulen : int32
        zulen : int32
        """
        _MiddleMan.f90wrap_middleman__mhdstartupsorted(xu=xu, yu=yu, zu=zu, \
            mhdposition_in=mhdposition_in, mhdb_in=mhdb_in, nx_split=nx_split, \
            ny_split=ny_split, nz_split=nz_split, mix=mix, max_bn=max_bn, miy=miy, \
            may=may, miz=miz, maz=maz, region_order_in=region_order_in, start_x=start_x, \
            end_x=end_x, start_y=start_y, end_y=end_y, start_z=start_z, end_z=end_z, \
            num_regions=num_regions, xulen=xulen, yulen=yulen, zulen=zulen)
    
    @staticmethod
    def gse2gswtsy15(date, position_gse, wind, gotso, hotso, glen, position_gsw, \
        interface_call=False):
        """
        gse2gswtsy15(date, position_gse, wind, gotso, hotso, glen, position_gsw)
        Defined at MiddleMan.f95 lines 1117-1147
        
        Parameters
        ----------
        date : float array
        position_gse : float array
        wind : float array
        gotso : float array
        hotso : float array
        glen : int32
        position_gsw : float array
        """
        _MiddleMan.f90wrap_middleman__gse2gswtsy15(date=date, position_gse=position_gse, \
            wind=wind, gotso=gotso, hotso=hotso, glen=glen, position_gsw=position_gsw)
    
    _dt_array_initialisers = []
    
    

middleman = Middleman()

