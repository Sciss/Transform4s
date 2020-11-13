package de.sciss.transform4s.utils

import scala.concurrent.{ExecutionContext, Future}

trait ConcurrencyUtilsPlatform {
  protected final var _executionContext: ExecutionContext =
    ExecutionContext.global

  /**
   * Waits for all threads to complete computation. Unsupported on Scala.js!
   */
  def waitForCompletion(futures: Array[Future[_]]): Unit =
    throw new UnsupportedOperationException
}