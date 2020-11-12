/* ***** BEGIN LICENSE BLOCK *****
 * JLargeArrays
 * Copyright (C) 2013 onward University of Warsaw, ICM
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ***** END LICENSE BLOCK ***** */

package de.sciss.transform4s.utils

import java.util.concurrent.Callable
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import java.util.concurrent.ThreadFactory
import java.util.concurrent.TimeUnit

import Math.max

/**
 * Concurrency utilities.
 * <p>
 *
 * @author Piotr Wendykier (p.wendykier@icm.edu.pl)
 */
object ConcurrencyUtils {
  /**
   * Thread pool.
   */
  private val DEFAULT_THREAD_POOL = Executors.newCachedThreadPool(new ConcurrencyUtils.CustomThreadFactory(new ConcurrencyUtils.CustomExceptionHandler))
  private var threadPool          = DEFAULT_THREAD_POOL
  private var nthreads            = getNumberOfProcessors
  private var concurrentThreshold = 100000L

  private class CustomExceptionHandler extends Thread.UncaughtExceptionHandler {
    override def uncaughtException(t: Thread, e: Throwable): Unit = {
      e.printStackTrace()
    }
  }

  private object CustomThreadFactory {
    private val DEFAULT_FACTORY = Executors.defaultThreadFactory
  }

  private class CustomThreadFactory(val handler: Thread.UncaughtExceptionHandler) extends ThreadFactory {
    override def newThread(r: Runnable): Thread = {
      val t = CustomThreadFactory.DEFAULT_FACTORY.newThread(r)
      t.setUncaughtExceptionHandler(handler)
      t
    }
  }

  /**
   * Returns the minimum length of array for which multiple threads are used.
   * <p>
   *
   * @return the minimum length of array for which multiple threads are used
   */
  def getConcurrentThreshold: Long = ConcurrencyUtils.concurrentThreshold

  /**
   * Sets the minimum length of an array for which multiple threads are used.
   * <p>
   *
   * @param concurrentThreshold minimum length of an array for which multiple threads are used
   */
  def setConcurrentThreshold(concurrentThreshold: Long): Unit = {
    ConcurrencyUtils.concurrentThreshold = max(1, concurrentThreshold)
  }

  /**
   * Returns the number of available processors.
   * <p>
   *
   * @return number of available processors
   */
  def getNumberOfProcessors: Int = Runtime.getRuntime.availableProcessors

  /**
   * Returns the current number of threads.
   * <p>
   *
   * @return the current number of threads.
   */
  def getNumberOfThreads: Int = ConcurrencyUtils.nthreads

  /**
   * Sets the number of threads.
   * <p>
   *
   * @param n new value of threads
   */
  def setNumberOfThreads(n: Int): Unit =
    ConcurrencyUtils.nthreads = n

  /**
   * Submits a value-returning task for execution and returns a Future
   * representing the pending results of the task.
   * <p>
   *
   * @tparam T  type
   * @param task task for execution
   *             <p>
   * @return handle to the task submitted for execution
   */
  def submit[T](task: Callable[T]): Future[T] = {
    if (ConcurrencyUtils.threadPool.isShutdown || ConcurrencyUtils.threadPool.isTerminated) {
      ConcurrencyUtils.threadPool = DEFAULT_THREAD_POOL
    }
    ConcurrencyUtils.threadPool.submit(task)
  }

  /**
   * Submits a Runnable task for execution and returns a Future representing that task.
   * <p>
   *
   * @param task task for execution
   *             <p>
   * @return handle to the task submitted for execution
   */
  def submit(task: Runnable): Future[_] = {
    if (ConcurrencyUtils.threadPool.isShutdown || ConcurrencyUtils.threadPool.isTerminated) ConcurrencyUtils.threadPool = DEFAULT_THREAD_POOL
    ConcurrencyUtils.threadPool.submit(task)
  }

  /**
   * Waits for all threads to complete computation.
   * <p>
   *
   * @param futures list of handles to the tasks
   *                <p>
   * @throws ExecutionException   if the computation threw an exception
   * @throws InterruptedException if the current thread was interrupted while waiting
   */
  @throws[InterruptedException]
  @throws[ExecutionException]
  def waitForCompletion(futures: Array[Future[_]]): Unit = {
    val size = futures.length
    for (j <- 0 until size) {
      futures(j).get
    }
  }

  /**
   * Sets the pool of threads.
   * <p>
   *
   * @param threadPool pool of threads
   */
  def setThreadPool(threadPool: ExecutorService): Unit = {
    ConcurrencyUtils.threadPool = threadPool
  }

  /**
   * Returns the pool of threads.
   * <p>
   *
   * @return pool of threads
   */
  def getThreadPool: ExecutorService = ConcurrencyUtils.threadPool

  /**
   * Shutdowns all submitted tasks.
   */
  def shutdownThreadPoolAndAwaitTermination(): Unit = {
    ConcurrencyUtils.threadPool.shutdown() // Disable new tasks from being submitted

    try // Wait a while for existing tasks to terminate
      if (!ConcurrencyUtils.threadPool.awaitTermination(60, TimeUnit.SECONDS)) {
        ConcurrencyUtils.threadPool.shutdownNow // Cancel currently executing tasks

        // Wait a while for tasks to respond to being cancelled
        if (!ConcurrencyUtils.threadPool.awaitTermination(60, TimeUnit.SECONDS)) System.err.println("Pool did not terminate")
      }
    catch {
      case _: InterruptedException =>
        // (Re-)Cancel if current thread also interrupted
        ConcurrencyUtils.threadPool.shutdownNow
        // Preserve interrupt status
        Thread.currentThread.interrupt()
    }
  }
}
