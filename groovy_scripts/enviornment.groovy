import org.jlab.utils.CLASResources;
import java.util.Optional;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.base.DetectorType;


String varname = CLASResources.getEnvironmentVariable("GEOMETRYDATABASEVARIATION");
// String varname = System.getenv("GEOMETRYDATABASEVARIATION");
String variationName = Optional.ofNullable(varname).orElse("default");
println variationName

ConstantProvider provider = GeometryFactory.getConstants(DetectorType.DC, 11, variationName);
// println provider.getDouble("/geometry/dc/region/dist2tgt",0)
println("the varnmae is "+varname);