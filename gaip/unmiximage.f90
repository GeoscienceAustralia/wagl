!-------------------------------------------------------
SUBROUTINE unmixImage(image, endMemberMatrix, numBands, numEndMembers, numRows, numCols, inNull, outNull, fractionsImage)
!-------------------------------------------------------
! Performs spectral unmixing of an image.
!
! For each pixel it finds the fractions that minimises
! || endMemberMatrix * fractions - pixelReflectance || ^ 2
!
! Uses the non-negative least squares (nnls) solver of Charles L. Lawson and
! Richard J. Hanson to do the optimization. See nnls.f90 for more details
! on the nnls algorithm.
!
! Tony Gill, January 2011.
!
! ARGUMENTS
! image (IN):              A 3D array containing the reflectance image (numBands, numRows, numCols).
! endMemberMatrix (IN):    A 2D array containing the reflectance of the endmembers (numBands, numEndmembers)
! numBands (IN):           The number of bands in the image.
! numEndMembers (IN):      The number of endmembers to be unmixed.
! numRows (IN):            The number of rows in the image.
! numCols (IN):            The number of columns in the image.
! inNull (IN):             The ignore value for the input image - if all values in the input pixel equal
!                          inNull then the pixel is not processed. The output fractions will be filled with outNull.
! outNull (IN):            The value to use to fill the output fractionsImage where the input pixel contains null data
!                          or the unmixing cannot be performed.
! fractionsImage (OUT):    A 3D array that on EXIT will contain an image.
!                          Each pixel in the image contains the endmember fractions plus
!                          an error term which is equal to the euclidean norm of the residual vector.
!                          The shape of the array is (numEndMembers+1,numRows,numCols).
!                          If nnls fails, the pixel will be filled with -10.0
!-------------------------------------------------------

   USE nonneg_leastsq
   IMPLICIT NONE

   ! Function arguments
   INTEGER, INTENT(IN)           :: numBands, numEndMembers, numRows, numCols
   DOUBLE PRECISION, INTENT(IN)  :: image(numBands,numRows,numCols), endMemberMatrix(numBands,numEndMembers)
   DOUBLE PRECISION, INTENT(IN)  :: inNull, outNull
   DOUBLE PRECISION, INTENT(OUT) :: fractionsImage(numEndMembers+1,numRows,numCols)

   ! Local variables
   INTEGER                       :: col, row, indx(numEndMembers), mode
   DOUBLE PRECISION              :: endMemberMatrixCopy(numBands,numEndMembers), pixelFractions(numEndMembers)
   DOUBLE PRECISION              :: pixelRefl(numBands), rnorm, w(numEndMembers)

   ! loop through the image.
   do col = 1, numCols
      do row = 1, numRows
        pixelRefl = image(:,row,col) ! pixel reflectance.

        if ( all(pixelRefl .EQ. inNull) ) then
           ! ignore the pixel - fill with out null.
           fractionsImage(:, row, col) = outNull
        else
           ! the pixel needs to be processed.
           endMemberMatrixCopy = endMemberMatrix ! nnls changes endMemberMatrix, so need to use a copy of it.
           CALL nnls(endMemberMatrixCopy, numBands, numEndMembers, pixelRefl, pixelFractions, rnorm, w, indx, mode)
           if (mode .GT. 1) then
              ! unmixing failed, fill the pixel.
              fractionsImage(:, row, col) = outNull
           else
              fractionsImage(:numEndMembers,row,col) = pixelFractions ! insert endmembers
              fractionsImage(numEndMembers+1, row, col) = rnorm ! insert rnorm value.
           endif
        endif
      end do
   end do
END SUBROUTINE unmixImage


!-------------------------------------------------------
SUBROUTINE unmixPixel(pixelRefl, endMemberMatrix, numBands, numEndMembers, inNull, outNull, pixelFractions)
!-------------------------------------------------------
! Performs spectral unmixing of a pixel.
!
! Finds the fractions that minimises
! || endMemberMatrix * fractions - pixelReflectance || ^ 2
!
! Uses the non-negative least squares (nnls) solver of Charles L. Lawson and
! Richard J. Hanson to do the optimization. See nnls.f90 for more details
! on the nnls algorithm.
!
! Tony Gill, January 2011.
!
! ARGUMENTS
! pixelRefl (IN):          A 1D array containing the reflectance (numBands).
! endMemberMatrix (IN):    A 2D array containing the reflectance of the endmembers (numBands, numEndmembers)
! numBands (IN):           The number of bands in the pixel.
! numEndMembers (IN):      The number of endmembers to be unmixed.
! inNull (IN):             The ignore value for the input image - if all values in the input pixel equal
!                          inNull then the pixel is not processed. The output fractions will be filled with outNull.
! outNull (IN):            The value to use to fill the output fractions array where the input pixel contains null data
!                          or the unmixing cannot be performed.
! pixelFractions (OUT):    A 1D array that on EXIT will contain the unmixed fractions.
!                          The array contains the endmember fractions plus
!                          an error term which is equal to the euclidean norm of the residual vector.
!                          The length of the array is (numEndMembers+1).
!-------------------------------------------------------

   USE nonneg_leastsq
   IMPLICIT NONE

   ! Function arguments
   INTEGER, INTENT(IN)           :: numBands, numEndMembers
   DOUBLE PRECISION, INTENT(IN)  :: pixelRefl(numBands), endMemberMatrix(numBands,numEndMembers)
   DOUBLE PRECISION, INTENT(IN)  :: inNull, outNull
   DOUBLE PRECISION, INTENT(OUT) :: pixelFractions(numEndMembers+1)

   ! Local variables
   INTEGER                       :: col, row, indx(numEndMembers), mode
   DOUBLE PRECISION              :: endMemberMatrixCopy(numBands,numEndMembers)
   DOUBLE PRECISION              :: pixelReflCopy(numBands), pixelFracs(numEndmembers), rnorm, w(numEndMembers)

   if ( all(pixelRefl .EQ. inNull) ) then
      ! ignore the pixel - fill with out null.
      pixelFractions(:) = outNull
   else
      ! the pixel needs to be processed.
      endMemberMatrixCopy = endMemberMatrix ! nnls changes endMemberMatrix, so need to use a copy of it.
      pixelReflCopy = pixelRefl
      CALL nnls(endMemberMatrixCopy, numBands, numEndMembers, pixelReflCopy, pixelFracs, rnorm, w, indx, mode)
      if (mode .GT. 1) then
         ! unmixing failed, fill the pixel.
         pixelFractions(:) = outNull
      else
         pixelFractions(:numEndMembers) = pixelFracs ! insert endmembers
         pixelFractions(numEndMembers+1) = rnorm ! insert rnorm value.
      endif
   endif
END SUBROUTINE unmixPixel
