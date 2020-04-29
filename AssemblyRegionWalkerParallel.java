package org.broadinstitute.hellbender.engine;

import com.google.common.util.concurrent.ThreadFactoryBuilder;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.filters.CountingReadFilter;
import org.broadinstitute.hellbender.transformers.ReadTransformer;
import org.broadinstitute.hellbender.utils.ResourcePool;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryUsage;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

public abstract class AssemblyRegionWalkerParallel extends AssemblyRegionWalker {
    @Argument(fullName = "useAutoAdjustment", doc = "auto adjustment calculate, default true")
    private boolean useAutoAdjustment = true;

    @Argument(fullName = "activeRegionThreads", doc = "use when useAutoAdjustment is false, active region threads", optional = true)
    private int activeRegionThreads = 1;

    @Argument(fullName = "callRegionThreads", doc = "use when useAutoAdjustment is false, call region threads", optional = true)
    private int callRegionThreads = 1;

    @Argument(fullName = "maxParallelThreads", doc = "use when useAutoAdjustment is true, " +
            "total parallel threads, used by auto parameter adjustment, default all available cpu cores", optional = true)
    private int maxParallelThreads = Runtime.getRuntime().availableProcessors()/2;

    @Argument(fullName = "init_threads", doc = "start contig numbers when enable useAutoAdjustment", optional = true)
    private int init_threads = maxParallelThreads >= 5 ? maxParallelThreads / 5 : 1;

    @Argument(fullName = "outputSchedulerStatus", doc = "whether output scheduler status", optional = true)
    private boolean outputSchedulerStatus = true;

    private ProgressMeterPar progressMeterPar;

    protected ResourcePool<InputsBundle> inputs = ResourcePool.withInitial(() -> {
        SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(readArguments.getReadValidationStringency());

        factory = factory.referenceSequence(referenceArguments.getReferencePath());

        if (bamIndexCachingShouldBeEnabled()) {
            factory = factory.enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES);
        }

        ReadsDataSource reads = new ReadsDataSource(readArguments.getReadPaths(),
                readArguments.getReadIndexPaths(), factory, cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer < 0 ? cloudPrefetchBuffer : cloudIndexPrefetchBuffer);

        ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());

        FeatureManager features = new FeatureManager(this, FeatureDataSource.DEFAULT_QUERY_LOOKAHEAD_BASES,
                cloudPrefetchBuffer, cloudIndexPrefetchBuffer, referenceArguments.getReferencePath());

        if (features.isEmpty()) {
            features = null;
        }

        return new InputsBundle(reads, reference, features);
    }).setCleaner((inputsBundle) -> {            //withInitial() will returns a SuppliedResourcePool object
        if (features != null)
            inputsBundle.features.close();
        inputsBundle.reads.close();
        inputsBundle.reference.close();
    });

    protected ResourcePool<FilterBundle> filters = ResourcePool.withInitial(() ->
            new FilterBundle(makeReadFilter(), createDownsampler(),
                    makePreReadFilterTransformer(), makePostReadFilterTransformer())
    );

    protected ResourcePool<AssemblyRegionEvaluator> enginePool = ResourcePool.withInitial(() -> newEngine())
            .setCleaner((engine) -> shutdownEngine(engine));

    protected final List<VariantContext> variantContexts = new ArrayList<>();

    private Scheduler scheduler;
    public static volatile boolean threadsPoolFinished = false;

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        progressMeterPar = new ProgressMeterPar(secondsBetweenProgressUpdates);
        progressMeterPar.setRecordLabel(getProgressMeterRecordLabel());
        scheduler = useAutoAdjustment ? new AutoAdjustmentCalculator() : new ManualAdjustmentCalculator();    //useAutoAdjustment=true
        scheduler.init();
        System.out.println("maxPArallelThreads: "+maxParallelThreads);
    }

    @Override
    public void traverse() {
        scheduler.start();
        scheduler.shutdown();
        threadsPoolFinished = true;
        System.out.println("Finish Successful, vcf size : " + variantContexts.size());
    }

    @Override
    public Object onTraversalSuccess() {
        progressMeterPar.stop();
        return null;
    }

    protected abstract AssemblyRegionEvaluator newEngine();

    protected abstract void shutdownEngine(AssemblyRegionEvaluator engine);

    @Override
    protected void processReadShard(MultiIntervalLocalReadShard shard, ReferenceDataSource reference, FeatureManager features) {
        final Iterator<AssemblyRegion> assemblyRegionIter = new AssemblyRegionIterator(shard, getHeaderForReads(), reference, features, assemblyRegionEvaluator(), minAssemblyRegionSize, maxAssemblyRegionSize, assemblyRegionPadding, activeProbThreshold, maxProbPropagationDistance, includeReadsWithDeletionsInIsActivePileups());

        while (assemblyRegionIter.hasNext()) {
            final AssemblyRegion assemblyRegion = assemblyRegionIter.next();

            // 预先过滤,避免创建过多线程
            if (!assemblyRegion.isActive() || assemblyRegion.size() == 0) {
                continue;
            }

            scheduler.addRegion(new CallRegionTask(assemblyRegion));
        }
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();
        inputs.shutdown();
    }

    private class ActiveRegionDetectTask implements Runnable, Comparable<ActiveRegionDetectTask> {
        private MultiIntervalLocalReadShard shard;

        ActiveRegionDetectTask(MultiIntervalLocalReadShard shard) {
            this.shard = shard;
        }

        @Override
        public void run() {
            try {
                // start process shard, register shard to progressMeter
                progressMeterPar.start(shard);
                // 多线程通过ThreadLocal获取各自独立的输入数据读取工具类,filter,以及engine
                InputsBundle inputsBundle = inputs.get();
                FilterBundle filterBundle = filters.get();
                shard.setReadsSource(inputsBundle.reads);

                //shard.setPreReadFilterTransformer(makePreReadFilterTransformer());
                shard.setPreReadFilterTransformer(filterBundle.preReadFilterTransformer);
                shard.setReadFilter(filterBundle.countedFilter);
                shard.setDownsampler(filterBundle.downsampler);
                shard.setPostReadFilterTransformer(filterBundle.postReadFilterTransformer);

                processReadShard(shard, inputsBundle.reference, inputsBundle.features);

                inputs.free();
                filters.free();
                scheduler.contigDone(shard);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        @Override
        public int compareTo(ActiveRegionDetectTask another) {
            int result = shard.getIntervals().size() - another.shard.getIntervals().size();
            if (result != 0)
                return result;
            else {
                int length1 = shard.getIntervals().stream().mapToInt((a) -> a.getEnd() - a.getStart()).sum();
                int length2 = another.shard.getIntervals().stream().mapToInt((a) -> a.getEnd() - a.getStart()).sum();
                return length1 - length2;
            }
        }
    }

    protected class CallRegionTask implements Runnable {
        AssemblyRegion region;

        public CallRegionTask(final AssemblyRegion region) {
            this.region = region;
        }

        @Override
        public void run() {
            try {
                InputsBundle inputsBundle = inputs.get();
                final ReferenceContext referenceContext = new ReferenceContext(inputsBundle.reference, region.getExtendedSpan());
                final FeatureContext featureContext = new FeatureContext(inputsBundle.features, region.getExtendedSpan());

                apply(region, referenceContext, featureContext);

                //free by pool
//                inputs.free();
//                enginePool.free();
                scheduler.regionDone(region);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    abstract class Scheduler {
        int contigsSize;
        int DEFAULT_SECONDS = (int) secondsBetweenProgressUpdates;
        PriorityQueue<ActiveRegionDetectTask> tasksQueue = new PriorityQueue<>(Comparator.reverseOrder());

        public Scheduler() {
            this.contigsSize = readShards.size();
            readShards.forEach((MultiIntervalLocalReadShard shard) -> tasksQueue.add(new ActiveRegionDetectTask(shard)));
        }

        abstract void init();

        abstract void start();

        abstract void addRegion(CallRegionTask task);

        void contigDone(MultiIntervalLocalReadShard shard) {
            progressMeterPar.contigDone(shard);
        }

        void regionDone(AssemblyRegion region) {
            if (progressMeterPar.update(region)) {
                adjustCalculator();
                if (outputSchedulerStatus)
                    System.out.print(getStatus());
            }
        }

        void adjustCalculator() {
        }

        // block method
        abstract void shutdown();

        abstract String getStatus();
    }

    class ManualAdjustmentCalculator extends Scheduler {
        private ThreadPoolExecutor activeRegionThreadsPool;
        private ThreadPoolExecutor callRegionThreadsPool;

        @Override
        public void init() {
            activeRegionThreadsPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(activeRegionThreads);
            callRegionThreadsPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(callRegionThreads,
                    new ThreadFactoryBuilder().setNameFormat("callRegion" + "-%d").build());
        }

        @Override
        public void start() {
            for (final ActiveRegionDetectTask task : tasksQueue) {
                // 这里的Filter和Transform应该保证线程安全 即使用函数式语言
                // 但是这里无法保证,所以后移到线程执行任务时进行创建,设定
                activeRegionThreadsPool.submit(task);
            }
        }

        @Override
        public void addRegion(CallRegionTask task) {
            callRegionThreadsPool.submit(task);
        }

        @Override
        public void shutdown() {
            activeRegionThreadsPool.shutdown();
            try {
                while (!activeRegionThreadsPool.awaitTermination(DEFAULT_SECONDS, TimeUnit.SECONDS)) {
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            callRegionThreadsPool.shutdown();
            try {
                while (!callRegionThreadsPool.awaitTermination(DEFAULT_SECONDS, TimeUnit.SECONDS)) {
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        @Override
        String getStatus() {
            return "[ThreadPool] active region active threads : " +
                    activeRegionThreadsPool.getActiveCount() +
                    " completed task count : " +
                    activeRegionThreadsPool.getCompletedTaskCount() +
                    " queue size : " +
                    activeRegionThreadsPool.getQueue().size() + "\n" +
                    "[ThreadPool] call region active threads : " +
                    callRegionThreadsPool.getActiveCount() +
                    " completed task count : " +
                    callRegionThreadsPool.getCompletedTaskCount() +
                    " queue size : " +
                    callRegionThreadsPool.getQueue().size() + "\n";
        }
    }

    class AutoAdjustmentCalculator extends Scheduler {           //class Scheduler contains a taskQueue which is actually a priority queue
        MemoryMXBean memoryMXBean = ManagementFactory.getMemoryMXBean();   //MemoryMXBean is the management interface for the memory system of java virtual machine
        private ThreadPoolExecutor threadsPool;

        AtomicInteger regionSize = new AtomicInteger();       //---these three AtomicInteger ?
        AtomicInteger regionComplete = new AtomicInteger();
        AtomicInteger contigsCompleted = new AtomicInteger();

        final Lock contigLock = new ReentrantLock();
        final Condition notFull = contigLock.newCondition();
        AtomicInteger INDEX = new AtomicInteger();           //---what's this ?

        Map<String, ContigStatus> processing = new ConcurrentHashMap<>();
        Map<String, ContigStatus> toBlock = new ConcurrentHashMap<>();        //mark which contig is blocked
        AtomicInteger blockedCount = new AtomicInteger();


        volatile long queueEmptyTime = System.currentTimeMillis();
        final Object queueEmptyTimeLock = new Object();
        volatile int queueSize;
        volatile int contigsProcessing = 0;

        class ContigStatus {
            String contig;
            int index;
            boolean canRunning;

            ContigStatus(String contig) {
                this.contig = contig;
                index = INDEX.getAndIncrement();
                canRunning = true;
            }
        }

        @Override
        public void init() {
            threadsPool = (ThreadPoolExecutor) Executors.newFixedThreadPool(maxParallelThreads,
                    new ThreadFactoryBuilder().setNameFormat("threadsPool" + "-%d").build());
        }

        @Override
        public void start() {
            for (int i = 0; i < init_threads; i++) {
                addContig(tasksQueue.remove());
            }
        }

        void addContig(ActiveRegionDetectTask task) {
            String name = task.shard.getIntervals().get(0).getContig();
            threadsPool.submit(task);                               //the only entrance to submit the task
            processing.put(name, new ContigStatus(name));
            contigsProcessing = processing.size();
            System.out.println("addContig: the number of active thread: " + threadsPool.getActiveCount());         //thread test
        }

        @Override
        void contigDone(MultiIntervalLocalReadShard shard) {
            super.contigDone(shard);
            String name = shard.getIntervals().get(0).getContig();
            processing.remove(name);
            contigsProcessing = processing.size();
            contigsCompleted.getAndIncrement();
            toBlock.remove(name);

            if (contigsProcessing - blockedCount.get() < init_threads && heapMemoryEnough()) {
                signalOrAddContig();
            }
            System.out.println("contigDone: the number of active thread: " + threadsPool.getActiveCount());         //thread test
        }

        @Override
        public void addRegion(CallRegionTask task) {
            regionSize.getAndIncrement();
            threadsPool.submit(task);
            if(regionSize.get()%1000==0)
                System.out.println("addRegion: the number of active thread: " + threadsPool.getActiveCount()+" regionSize:"+regionSize.get());         //thread test
            ContigStatus status;
            if ((status = toBlock.get(task.region.getContig())) != null) {
                System.out.println("[Auto adjustment] block " + task.region.getContig() + " id :" + status.index);
                contigLock.lock();
                try {
                    blockedCount.getAndIncrement();
                    while (queueSize > maxParallelThreads) {
                        notFull.await();        //---await() method will block the current thread
                    }
                } catch (InterruptedException e) {
                    e.printStackTrace();
                } finally {
                    contigLock.unlock();
                }

                toBlock.remove(task.region.getContig());
                status.canRunning = true;
                blockedCount.getAndDecrement();
                System.out.println("[Auto adjustment] weak up " + task.region.getContig() + " id :" + status.index);
            }
        }

        @Override
        void regionDone(AssemblyRegion region) {
            super.regionDone(region);
            regionComplete.getAndIncrement();
        }

        @Override
        void adjustCalculator() {
            queueSize = threadsPool.getQueue().size();

            if (shouldAddContig()) {
                signalOrAddContig();
                return;
            }

            if (shouldMinusContig()) {
                markContig();
            }
        }

        private boolean shouldMinusContig() {
            return !heapMemoryEnough() && contigsProcessing - blockedCount.get() > init_threads;
        }


        boolean heapMemoryEnough() {
            MemoryUsage memoryUsage = memoryMXBean.getHeapMemoryUsage();
            return (memoryUsage.getUsed() / (double) memoryUsage.getMax()) < 0.9;
        }

        private void markContig() {
            int max = -1;
            String name = null;
            ContigStatus status = null;
            for (Map.Entry<String, ContigStatus> entry : processing.entrySet()) {    //---record the processing contigs using map
                if (entry.getValue().canRunning && entry.getValue().index > max) {
                    max = entry.getValue().index;
                    name = entry.getKey();
                    status = entry.getValue();
                }
            }

            if (name != null && toBlock.get(name) == null) {
                status.canRunning = false;
                toBlock.put(name, status);
                System.out.println("[Auto adjustment] markContig : " + name);
            }
        }
                                                          //this method is called everytime a region is done
        private boolean shouldAddContig() {               //determine whether it's necessary to add a thread for a new readShard
            if (contigsCompleted.get() + contigsProcessing - blockedCount.get() == contigsSize)    //---there is no contigs which haven't been dealed with
                return false;
            if (threadsPool.getActiveCount() < maxParallelThreads)
                return true;

            if (queueSize == 0) {                           //callRegion task queue is empty    //---???
                synchronized (queueEmptyTimeLock) {
                    long current = System.currentTimeMillis();
                    if (current - queueEmptyTime > 2000) {
                        queueEmptyTime = current;
                        return true;
                    }
                }
            }
            return false;
        }

        private void signalOrAddContig() {      //if there is any contig blocked, then signal to awake it, if not ,then add a contig
            if (blockedCount.get() > 0) {
                contigLock.lock();
                try {
                    notFull.signal();
                } finally {
                    contigLock.unlock();
                }
                return;
            }

            if (tasksQueue.size() > 0) {
                contigLock.lock();
                if (tasksQueue.size() > 0) {
                    try {
                        addContig(tasksQueue.remove());
                    } finally {
                        contigLock.unlock();
                    }
                }
            }
        }

        @Override
        public void shutdown() {
            while (contigsCompleted.get() < contigsSize) {       //contigsSize=86  the size of readShard
                try {
                    Thread.sleep(10 * 1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            threadsPool.shutdown();
            try {
                while (!threadsPool.awaitTermination(DEFAULT_SECONDS, TimeUnit.SECONDS)) {
                }
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }

        @Override
        String getStatus() {
            return "contig : " + contigsCompleted + " / " + contigsSize + " processing : " + contigsProcessing + " blocked : " + toBlock.size() + "\n"
                    + "region : " + regionComplete + " / " + regionSize + " queue size : " + queueSize + "\n";
        }
    }

    public static class FilterBundle {
        CountingReadFilter countedFilter;
        ReadsDownsampler downsampler;
        ReadTransformer preReadFilterTransformer;
        ReadTransformer postReadFilterTransformer;

        public FilterBundle(CountingReadFilter countedFilter, ReadsDownsampler downsampler, ReadTransformer preReadFilterTransformer, ReadTransformer postReadFilterTransformer) {
            this.countedFilter = countedFilter;
            this.downsampler = downsampler;
            this.preReadFilterTransformer = preReadFilterTransformer;
            this.postReadFilterTransformer = postReadFilterTransformer;
        }
    }

    public static class InputsBundle {
        public ReadsDataSource reads;
        public ReferenceDataSource reference;
        public FeatureManager features;

        InputsBundle(ReadsDataSource reads, ReferenceDataSource reference, FeatureManager features) {
            this.reads = reads;
            this.reference = reference;
            this.features = features;
        }
    }
}
