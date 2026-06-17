module mpas_output_timelevels

    use mpas_kind_types
    use mpas_derived_types
    use mpas_timekeeping
    use mpas_log
    use mpas_string_utils, only : mpas_split_string


    public :: parse_all_timelevels, &
              get_output_interval_from_timelevels, &
              get_next_timelevel_start, &
              update_variable_output_alarm

    private


    contains


    !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !
    !  parse_all_timelevels
    !
    !> \brief Parse the entire output_timelevels string into a timelevel_spec structure
    !> \author Guoqing Ge
    !> \date March 2026
    !> \details
    !>  Parses the complete output_timelevels string (e.g., "0-3-15m 4-72 75-168-3")
    !>  into a pre-parsed structure for efficient runtime lookup. Uses mpas_split_string
    !>  to tokenize the string by spaces. Reports parse errors via mpas_log_write.
    !
    !-----------------------------------------------------------------------
    subroutine parse_all_timelevels(timelevels_str, timelevel_spec, stream_name, ierr)!{{{

        implicit none

        character(len=*), intent(in) :: timelevels_str
        type(MPAS_timelevel_spec_type), intent(out) :: timelevel_spec
        character(len=*), intent(in) :: stream_name
        integer, intent(out) :: ierr

        character(len=StrKIND), pointer, dimension(:) :: segments
        character(len=StrKIND) :: segment
        integer :: i, n_segments, local_ierr
        real(kind=RKIND) :: start_hour, end_hour, interval_mins

        ierr = 0
        timelevel_spec % n_segments = 0
        timelevel_spec % is_parsed = .false.

        ! Guard against empty string
        if (len_trim(timelevels_str) == 0) then
            call mpas_log_write('Stream '''//trim(stream_name)//''': output_timelevels string is empty', MPAS_LOG_ERR)
            ierr = 1
            return
        end if

        ! Split the string by spaces (trim leading/trailing whitespaces first)
        nullify(segments)
        call mpas_split_string(trim(adjustl(timelevels_str)), ' ', segments)
        n_segments = size(segments)

        ! Check for too many segments
        if (n_segments > MAX_TIMELEVEL_SEGMENTS) then
            call mpas_log_write('Stream '''//trim(stream_name)//''': too many timelevel segments (max $i)', &
                               MPAS_LOG_ERR, intArgs=(/MAX_TIMELEVEL_SEGMENTS/))
            deallocate(segments)
            ierr = 1
            return
        end if

        ! Parse each segment
        do i = 1, n_segments
            segment = trim(adjustl(segments(i)))

            ! Skip empty segments (from consecutive spaces)
            if (len_trim(segment) == 0) cycle

            call parse_output_timelevel_spec(segment, start_hour, end_hour, interval_mins, local_ierr)

            if (local_ierr /= 0) then
                call mpas_log_write('Stream '''//trim(stream_name)//''': Invalid timelevel format. ' // &
                                   'Expected ''start-end-interval'' but got '''//trim(segment)//'''', MPAS_LOG_ERR)
                deallocate(segments)
                ierr = 1
                return
            end if

            ! Store parsed values
            timelevel_spec % n_segments = timelevel_spec % n_segments + 1
            timelevel_spec % start_hour(timelevel_spec % n_segments) = start_hour
            timelevel_spec % end_hour(timelevel_spec % n_segments) = end_hour
            timelevel_spec % interval_minutes(timelevel_spec % n_segments) = interval_mins
        end do

        deallocate(segments)

        ! Mark as successfully parsed if we have at least one segment
        if (timelevel_spec % n_segments > 0) then
            timelevel_spec % is_parsed = .true.
        else
            call mpas_log_write('Stream '''//trim(stream_name)//''': no valid timelevel segments found', MPAS_LOG_ERR)
            ierr = 1
        end if

    end subroutine parse_all_timelevels!}}}


    !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !
    !  parse_time_string
    !
    !> \brief Parse a time string into total minutes ( allow 0.5 minutes, i.e. 30 seconds)
    !> \author Guoqing Ge
    !> \date February 2026
    !> \details
    !>  Parses time strings like "1h30m", "45m", "6", "90s", "1d" into total minutes
    !>  Supported units: d/D (days), h (hours), m (minutes), s (seconds)
    !>  Plain integers without units are interpreted as hours
    !
    !-----------------------------------------------------------------------
    subroutine parse_time_string(time_str, total_minutes, ierr)!{{{

        implicit none

        character(len=*), intent(in) :: time_str
        real(kind=RKIND), intent(out) :: total_minutes
        integer, intent(out) :: ierr

        character(len=StrKIND) :: str_local
        integer :: i, len_str, num_start, read_err
        real(kind=RKIND) :: value
        character :: ch

        ierr = 0
        total_minutes = 0.0_RKIND
        str_local = trim(adjustl(time_str))
        len_str = len_trim(str_local)

        if (len_str == 0) then
            ierr = 1
            return
        end if

        ! Check if last character is a digit (plain integer = hours)
        ch = str_local(len_str:len_str)
        if (ch >= '0' .and. ch <= '9') then
            read(str_local, *, iostat=read_err) value
            if (read_err /= 0) then
                ierr = 1
                return
            end if
            total_minutes = value * 60.0_RKIND
            return
        end if

        ! Parse duration with units (e.g., "1h30m", "45m", "90s", "1d" or "1D")
        num_start = 1
        i = 1
        do while (i <= len_str)
            ch = str_local(i:i)
            ! Check for unit characters: d/D (days), h (hours), m (minutes), s (seconds)
            if (ch == 'h' .or. ch == 'm' .or. ch == 's' .or. ch == 'd' .or. ch == 'D') then
                
                if (i == num_start) then
                    ierr = 1  ! No number before unit
                    return
                end if
                
                read(str_local(num_start:i-1), *, iostat=read_err) value
                if (read_err /= 0) then
                    ierr = 1
                    return
                end if

                select case (ch)
                    case ('h')
                        total_minutes = total_minutes + value * 60.0_RKIND
                    case ('m')
                        total_minutes = total_minutes + value
                    case ('s')
                        total_minutes = total_minutes + value / 60.0_RKIND
                    case ('d', 'D')
                        total_minutes = total_minutes + value * 24.0_RKIND * 60.0_RKIND
                end select

                num_start = i + 1
            end if
            i = i + 1
        end do

        ! Check for trailing digits without unit (error)
        if (num_start <= len_str) then
            ierr = 1
            return
        end if

    end subroutine parse_time_string!}}}


    !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !
    !  parse_output_timelevel_spec
    !
    !> \brief Parse one time level specification segment
    !> \author Guoqing Ge
    !> \date February 2026
    !> \details
    !>  Parses a single segment using format: start, start-stop, or start-stop-step
    !>  Time strings can be integers (hours) or duration format (e.g., "1h30m", "45m")
    !>  Examples: "6", "0-3", "0-1h-15m", "1h30m-2h-15m"
    !
    !-----------------------------------------------------------------------
    subroutine parse_output_timelevel_spec(spec, start_hour, end_hour, interval_minutes, ierr)!{{{

        implicit none

        character(len=*), intent(in) :: spec
        real(kind=RKIND), intent(out) :: start_hour, end_hour
        real(kind=RKIND), intent(out) :: interval_minutes
        integer, intent(out) :: ierr

        character(len=StrKIND), pointer, dimension(:) :: parts
        integer :: i, n_parts, local_ierr
        real(kind=RKIND) :: start_minutes, end_minutes, step_minutes

        ierr = 0

        ! Guard against empty string
        if (len_trim(spec) == 0) then
            ierr = 1
            return
        end if

        ! Split by dash delimiter
        nullify(parts)
        call mpas_split_string(trim(adjustl(spec)), '-', parts)
        n_parts = size(parts)

        ! Validate parts array - check for empty parts (from consecutive dashes or leading/trailing dashes)
        do i = 1, n_parts
            if (len_trim(parts(i)) == 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
        end do

        ! Parse based on number of parts
        if (n_parts == 1) then
            ! Format: single time string (output at that time only)
            call parse_time_string(parts(1), start_minutes, local_ierr)
            if (local_ierr /= 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
            start_hour = start_minutes / 60.0_RKIND
            end_hour = start_hour
            interval_minutes = 60.0_RKIND  ! Default, but won't matter for single time

        else if (n_parts == 2) then
            ! Format: start-stop (interval defaults to 1 hour)
            call parse_time_string(parts(1), start_minutes, local_ierr)
            if (local_ierr /= 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
            
            call parse_time_string(parts(2), end_minutes, local_ierr)
            if (local_ierr /= 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
            
            start_hour = start_minutes / 60.0_RKIND
            end_hour = end_minutes / 60.0_RKIND
            interval_minutes = 60.0_RKIND

        else if (n_parts == 3) then
            ! Format: start-stop-step
            call parse_time_string(parts(1), start_minutes, local_ierr)
            if (local_ierr /= 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
            
            call parse_time_string(parts(2), end_minutes, local_ierr)
            if (local_ierr /= 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
            
            call parse_time_string(parts(3), step_minutes, local_ierr)
            if (local_ierr /= 0) then
                ierr = 1
                deallocate(parts)
                return
            end if
            
            start_hour = start_minutes / 60.0_RKIND
            end_hour = end_minutes / 60.0_RKIND
            interval_minutes = step_minutes

        else
            ! Invalid number of parts (0 or > 3 dashes)
            ierr = 1
            deallocate(parts)
            return
        end if

        deallocate(parts)

        ! Check for setting mistake "start-step-stop". The correct one should be "start-stop-step"
        if ((n_parts == 3) .and. (interval_minutes > end_hour * 60.0_RKIND + 1.0e-4_RKIND)) then
            ierr = 1
            call mpas_log_write('ERROR: output_timelevels segment '''//trim(spec)//''': ' // &
                               'expected format is start-stop-step, but appears to be start-step-stop', MPAS_LOG_ERR)
        end if

    end subroutine parse_output_timelevel_spec!}}}


    !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !
    !  get_output_interval_from_timelevels
    !
    !> \brief Calculate output interval based on current forecast hour and pre-parsed timelevels
    !> \author Guoqing Ge
    !> \date March 2026
    !> \details
    !>  Given a pre-parsed timelevel_spec and current forecast hour, determines the
    !>  appropriate output interval in minutes.
    !
    !-----------------------------------------------------------------------
    subroutine get_output_interval_from_timelevels(timelevel_spec, forecast_hour, interval_minutes, ierr, next_output_hour)!{{{

        implicit none

        type(MPAS_timelevel_spec_type), intent(in) :: timelevel_spec
        real(kind=RKIND), intent(in) :: forecast_hour
        real(kind=RKIND), intent(out) :: interval_minutes
        integer, intent(out) :: ierr
        real(kind=RKIND), intent(out), optional :: next_output_hour

        integer :: i
        real(kind=RKIND) :: start_hour, end_hour, seg_interval
        real(kind=RKIND) :: next_hour
        logical :: found, is_single_time

        ierr = 0
        interval_minutes = 60.0_RKIND  ! Default to 1 hour
        if (present(next_output_hour)) next_output_hour = -1.0_RKIND
        found = .false.
        is_single_time = .false.
        next_hour = huge(1.0_RKIND)

        ! Guard against unparsed spec
        if (.not. timelevel_spec % is_parsed .or. timelevel_spec % n_segments == 0) then
            ierr = 1
            return
        end if

        ! First pass: find if current forecast hour matches any segment
        do i = 1, timelevel_spec % n_segments
            start_hour = timelevel_spec % start_hour(i)
            end_hour = timelevel_spec % end_hour(i)
            seg_interval = timelevel_spec % interval_minutes(i)

            ! Check if current forecast hour falls in this range
            if (forecast_hour >= start_hour .and. forecast_hour <= end_hour) then
                found = .true.
                ! Check if this is a single time point (not a range)
                if (abs(start_hour - end_hour) < 1.0e-6_RKIND) then
                    is_single_time = .true.
                else
                    ! It's a range - use the segment's interval
                    interval_minutes = seg_interval
                    return
                end if
            end if
        end do

        ! If we matched a single time, find the next output time
        if (found .and. is_single_time) then
            do i = 1, timelevel_spec % n_segments
                start_hour = timelevel_spec % start_hour(i)
                end_hour = timelevel_spec % end_hour(i)

                ! For single times, check if > current
                if (abs(start_hour - end_hour) < 1.0e-6_RKIND) then
                    if (start_hour > forecast_hour + 1.0e-6_RKIND .and. start_hour < next_hour) then
                        next_hour = start_hour
                    end if
                else
                    ! For ranges, check if start is after current
                    if (start_hour > forecast_hour + 1.0e-6_RKIND .and. start_hour < next_hour) then
                        next_hour = start_hour
                    end if
                end if
            end do

            ! Calculate interval to next time
            if (next_hour < huge(1.0_RKIND)) then
                interval_minutes = (next_hour - forecast_hour) * 60.0_RKIND
                if (present(next_output_hour)) next_output_hour = next_hour
            else
                ! No more times after this - interval_minutes = 0 signals last output
                interval_minutes = 0.0_RKIND
                if (present(next_output_hour)) next_output_hour = -1.0_RKIND
            end if
            return
        end if

        ! If no matching range found, signal that no more output is needed
        if (.not. found) then
            interval_minutes = 0.0_RKIND
            ierr = MPAS_STREAM_MGR_ERROR  ! Signal no more output
        end if

    end subroutine get_output_interval_from_timelevels!}}}


    !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !
    !  get_next_timelevel_start
    !
    !> \brief Find the start hour of the next timelevel range after current forecast hour
    !> \author Guoqing Ge
    !> \date February 2026
    !> \details
    !>  Given a time levels string and current forecast hour, finds the start hour
    !>  of the next timelevel range (if any). Returns -1 if no future ranges exist.
    !
    !-----------------------------------------------------------------------
    subroutine get_next_timelevel_start(timelevel_spec, forecast_hour, next_start_hour, next_interval_minutes, ierr)!{{{

        implicit none

        type(MPAS_timelevel_spec_type), intent(in) :: timelevel_spec
        real(kind=RKIND), intent(in) :: forecast_hour
        real(kind=RKIND), intent(out) :: next_start_hour
        real(kind=RKIND), intent(out) :: next_interval_minutes
        integer, intent(out) :: ierr

        integer :: i
        real(kind=RKIND) :: start_hour, best_start, best_interval
        logical :: found

        ierr = 0
        next_start_hour = -1.0_RKIND
        next_interval_minutes = 60.0_RKIND
        best_start = huge(1.0_RKIND)
        best_interval = 60.0_RKIND
        found = .false.

        ! Guard against unparsed spec
        if (.not. timelevel_spec % is_parsed .or. timelevel_spec % n_segments == 0) then
            ierr = 1
            return
        end if

        ! Find the segment with smallest start_hour > forecast_hour
        do i = 1, timelevel_spec % n_segments
            start_hour = timelevel_spec % start_hour(i)

            ! Check if this range starts after current forecast hour
            if (start_hour > forecast_hour .and. start_hour < best_start) then
                best_start = start_hour
                best_interval = timelevel_spec % interval_minutes(i)
                found = .true.
            end if
        end do

        if (found) then
            next_start_hour = best_start
            next_interval_minutes = best_interval
        else
            ierr = MPAS_STREAM_MGR_ERROR  ! No future ranges
        end if

    end subroutine get_next_timelevel_start!}}}


    !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    !
    !  update_variable_output_alarm
    !
    !> \brief Update variable output alarm with new interval based on forecast hour
    !> \author Guoqing Ge
    !> \date March 2026
    !> \details
    !>  This routine updates a variable output alarm after it has rung by:
    !>  1. Checking if the stream uses output_timelevels (via timelevel_spec)
    !>  2. Calculating current forecast hour (current time - start time)
    !>  3. Getting the next interval from get_output_interval_from_timelevels
    !>  4. Removing the old alarm and creating a new one with the updated interval
    !
    !-----------------------------------------------------------------------
    subroutine update_variable_output_alarm(manager, stream, alarm_name, ierr)!{{{

        implicit none

        type (MPAS_streamManager_type), intent(inout) :: manager
        type (MPAS_stream_list_type), pointer :: stream
        character(len=*), intent(in) :: alarm_name
        integer, intent(out), optional :: ierr

        type (MPAS_Time_type) :: start_time, current_time, alarmTime_local
        type (MPAS_TimeInterval_type) :: time_diff, alarmInterval_local
        real (kind=RKIND) :: forecast_hour, next_start_hour, next_output_hour
        real (kind=RKIND) :: interval_minutes, next_interval_minutes
        integer :: local_ierr, ierr_tmp
        integer (kind=I8KIND) :: seconds_diff

        local_ierr = MPAS_STREAM_MGR_NOERR

        ! Check if stream uses output_timelevels (via pre-parsed timelevel_spec)
        if (.not. stream % timelevel_spec % is_parsed) then
            ! Not a variable output stream, nothing to do
            if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
            return
        end if

        ! Calculate current forecast hour
        start_time = mpas_get_clock_time(manager % streamClock, MPAS_START_TIME, ierr=ierr_tmp)
        current_time = mpas_get_clock_time(manager % streamClock, MPAS_NOW, ierr=ierr_tmp)
        time_diff = current_time - start_time
        call mpas_get_timeInterval(time_diff, S_i8=seconds_diff, ierr=ierr_tmp)
        forecast_hour = real(seconds_diff, RKIND) / 3600.0_RKIND

        ! Get next interval using pre-parsed timelevel_spec
        ! Also get the absolute next output hour to avoid floating-point drift
        call get_output_interval_from_timelevels(stream % timelevel_spec, forecast_hour, interval_minutes, ierr_tmp, next_output_hour)
        if (ierr_tmp /= MPAS_STREAM_MGR_NOERR .or. interval_minutes <= 0.0_RKIND) then
            ! Current hour not in any range - check for future range
            call get_next_timelevel_start(stream % timelevel_spec, forecast_hour, next_start_hour, next_interval_minutes, ierr_tmp)
            call mpas_remove_clock_alarm(manager % streamClock, alarm_name, ierr=ierr_tmp)
            
            if (next_start_hour > 0.0_RKIND) then
                ! Schedule recurring alarm for when the next range starts
                call mpas_set_timeInterval(time_diff, dt=next_start_hour * 3600.0_RKIND, ierr=ierr_tmp)
                alarmTime_local = start_time + time_diff
                call mpas_set_timeInterval(alarmInterval_local, dt=next_interval_minutes * 60.0_RKIND, ierr=ierr_tmp)
                call mpas_add_clock_alarm(manager % streamClock, alarm_name, alarmTime_local, &
                                         alarmTimeInterval=alarmInterval_local, ierr=ierr_tmp)
            end if
            if (present(ierr)) ierr = MPAS_STREAM_MGR_NOERR
            return
        end if

        ! Convert minutes to time interval (using dt in seconds for sub-minute precision)
        call mpas_set_timeInterval(alarmInterval_local, dt=interval_minutes * 60.0_RKIND, ierr=ierr_tmp)
        local_ierr = ior(local_ierr, ierr_tmp)

        ! Remove old alarm
        call mpas_remove_clock_alarm(manager % streamClock, alarm_name, ierr=ierr_tmp)
        local_ierr = ior(local_ierr, ierr_tmp)

        ! Add new alarm with updated interval
        ! Use ABSOLUTE time from start_time to avoid floating-point drift
        ! This ensures discrete times like "0m 15m 1h3m 2 3" hit exactly
        if (next_output_hour > 0.0_RKIND) then
            ! Compute alarm time absolutely: start_time + next_output_hour
            call mpas_set_timeInterval(time_diff, dt=next_output_hour * 3600.0_RKIND, ierr=ierr_tmp)
            alarmTime_local = start_time + time_diff
        else
            ! Fallback to relative calculation (for ranges)
            alarmTime_local = current_time + alarmInterval_local
        end if
        call mpas_add_clock_alarm(manager % streamClock, alarm_name, alarmTime_local, &
                                 alarmTimeInterval=alarmInterval_local, ierr=ierr_tmp)
        local_ierr = ior(local_ierr, ierr_tmp)

        if (present(ierr)) ierr = local_ierr

    end subroutine update_variable_output_alarm!}}}


end module mpas_output_timelevels
