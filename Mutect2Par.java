package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.annotator.*;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;
import org.broadinstitute.hellbender.utils.variant.writers.SomaticGVCFWriter;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

/**
 * 多线程并发Mutect2
 */
@CommandLineProgramProperties(
        summary = "Call somatic SNVs and indels via local assembly of haplotypes",
        oneLineSummary = "Call somatic SNVs and indels via local assembly of haplotypes",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class Mutect2Par extends AssemblyRegionWalkerParallel {

    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public File outputVCF;

    private VariantContextWriter vcfWriter;

    private Mutect2Engine m2Engine;

    @Override
    public AssemblyRegionEvaluator newEngine(){
        VariantAnnotatorEngine annotatorEngine = new VariantAnnotatorEngine(makeVariantAnnotations(), null, Collections.emptyList(), false, false);
        return new Mutect2Engine(MTAC, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), referenceArguments.getReferenceFileName(), annotatorEngine);

    }

    @Override
    public void shutdownEngine(AssemblyRegionEvaluator engine){
        ((Mutect2Engine)engine).shutdown();
    }

    @Override
    protected int defaultMinAssemblyRegionSize() {
        return 50;
    }

    @Override
    protected int defaultMaxAssemblyRegionSize() {
        return 300;
    }

    @Override
    protected int defaultAssemblyRegionPadding() {
        return 100;
    }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() {
        return 50;
    }

    @Override
    protected double defaultActiveProbThreshold() {
        return 0.002;
    }

    @Override
    protected int defaultMaxProbPropagationDistance() {
        return 50;
    }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() {
        return true;
    }

    @Override
    public boolean useVariantAnnotations() { return true;}

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    public ReadTransformer makePostReadFilterTransformer() {
        return super.makePostReadFilterTransformer().andThen(Mutect2Engine.makeStandardMutect2PostFilterReadTransformer(referenceArguments.getReferencePath(), !MTAC.dontClipITRArtifacts));
    }

    @Override
    public List<Class<? extends Annotation>> getDefaultVariantAnnotationGroups() {
        return Mutect2Engine.getStandardMutect2AnnotationGroups();
    }

    @Override
    protected ReadsDownsampler createDownsampler() {
        return maxReadsPerAlignmentStart > 0 ?
                new MutectDownsampler(maxReadsPerAlignmentStart, MTAC.maxSuspiciousReadsPerAlignmentStart, MTAC.downsamplingStride) : null;
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() {
        return enginePool.get();
    }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        m2Engine = (Mutect2Engine)enginePool.get();
        if (m2Engine.emitReferenceConfidence()) {
            logger.warn("Note that the Mutect2 reference confidence mode is in BETA -- the likelihoods model and output format are subject to change in subsequent versions.");
            if ( MTAC.emitReferenceConfidence == ReferenceConfidenceMode.GVCF ) {
                try {
                    vcfWriter = new SomaticGVCFWriter(vcfWriter, new ArrayList<Number>(MTAC.GVCFGQBands));
                } catch ( IllegalArgumentException e ) {
                    throw new CommandLineException.BadArgumentValue("GQBands", "are malformed: " + e.getMessage());
                }
            }
        }


        vcfWriter = createVCFWriter(outputVCF);
      //  vcfWriter = createVCFWriter(outputVCF);
        m2Engine.writeHeader(vcfWriter, getDefaultToolVCFHeaderLines());
    }

    @Override
    public Collection<Annotation> makeVariantAnnotations(){
        final Collection<Annotation> annotations = super.makeVariantAnnotations();

        if (MTAC.mitochondria) {
            annotations.add(new OriginalAlignment());
        }
        return annotations;
    }

    @Override                                       //apply method is for AssemblyRegion
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        List<VariantContext> result = ((Mutect2Engine)enginePool.get()).callRegion(region, referenceContext, featureContext);
        synchronized (variantContexts) {                //---pay attention here !
            variantContexts.addAll(result);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        writeVcf();
        return "SUCCESS";
    }

    private void writeVcf() {
        variantContexts.sort(m2Engine.makeVCFHeader(Collections.emptySet()).getVCFRecordComparator());
        variantContexts.forEach(vcfWriter::add);
    }

    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
        enginePool.shutdown();
    }
}
