function fasores()
    % Clear workspace and figures
    clear all;
    close all;
    clc;

    % Display welcome message
    fprintf('Circuit Impedance Calculator\n');
    fprintf('---------------------------\n');

    % Get user input with validation
    [VS, VR, VZL, R] = get_validated_input();

    % Calculate intersection points of circles and phasors
    [x_intersect, y_intersect, angle_VR, angle_VZL] = calculate_fasors(VS, VR, VZL);

    % Plot phasor diagram
    plot_phasor_diagram(VS, VR, VZL, x_intersect, y_intersect);

    % Calculate and display impedance
    calculate_and_display_impedance(VS, VR, VZL, R, angle_VZL);
end

% Function to get and validate user input
function [VS, VR, VZL, R] = get_validated_input()
    valid_input = false;

    while ~valid_input
        fprintf('\nEnter circuit parameters:\n');
        VS = input('Source voltage (VS) in Volts: ');
        VR = input('Resistor voltage (VR) in Volts: ');
        VZL = input('Load impedance voltage (VZL) in Volts: ');
        R = input('Resistor value (R) in Ohms: ');

        % Validate input: VS must be less than or equal to VR + VZL
        if VS > (VR + VZL)
            fprintf('\nError: VS cannot be greater than VR + VZL\n');
            fprintf('Please enter values again.\n');
        else
            valid_input = true;
        end
    end
end

% Function to calculate intersection points of circles
function [x_intersect, y_intersect, angle_VR, angle_VZL] = calculate_fasors(VS, VR, VZL)
    % Solve for intersection of two circles:
    % Circle 1: x^2 + y^2 = VR^2 (centered at origin)
    % Circle 2: (x-VS)^2 + y^2 = VZL^2 (centered at (VS,0))

    % Using algebraic solution
    d = VS; % distance between circle centers

    % Check if circles intersect
    if d > VR + VZL || d < abs(VR - VZL)
        error('Circles do not intersect with these parameters');
    end

    % Calculate intersection points using the formula
    a = (VR^2 - VZL^2 + d^2) / (2 * d);
    h = sqrt(VR^2 - a^2);

    % Calculate intersection points
    x2 = a; % Since we're aligning VS along x-axis
    y2_pos = h; % positive y solution
    y2_neg = -h; % negative y solution

    % For the positive quadrant (as requested)
    x_intersect = x2;
    y_intersect = y2_pos;

    % Calculate angles
    angle_VR = atan2(y_intersect, x_intersect);
    angle_VZL = atan2(-y_intersect, VS - x_intersect);

    fprintf('Intersection point: (%.4f, %.4f)\n', x_intersect, y_intersect);
    fprintf('VR angle: %.4f degrees\n', rad2deg(angle_VR));
    fprintf('VZL angle: %.4f degrees\n', rad2deg(angle_VZL));
end

% Function to plot phasor diagram
function plot_phasor_diagram(VS, VR, VZL, x_intersect, y_intersect)
    figure;
    hold on;
    grid on;
    axis equal;

    % Plot the circles
    theta = linspace(0, 2*pi, 100);

    % Circle 1 (VR)
    x_circle1 = VR * cos(theta);
    y_circle1 = VR * sin(theta);
    plot(x_circle1, y_circle1, 'b--');

    % Circle 2 (VZL)
    x_circle2 = VS + VZL * cos(theta);
    y_circle2 = VZL * sin(theta);
    plot(x_circle2, y_circle2, 'r--');

    % Plot the fasors
    quiver(0, 0, x_intersect, y_intersect, 0, 'b', 'LineWidth', 2); % VR fasor
    quiver(VS, 0, x_intersect-VS, y_intersect, 0, 'r', 'LineWidth', 2); % VZL fasor
    quiver(0, 0, VS, 0, 0, 'k', 'LineWidth', 2); % VS fasor

    % Plot the current vector (perpendicular to voltage)
    current_length = 0.5; % for visualization
    quiver(x_intersect, y_intersect, -y_intersect*current_length, x_intersect*current_length, 0, 'g', 'LineWidth', 2);

    % Mark the intersection point
    plot(x_intersect, y_intersect, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    % Labels
    text(x_intersect/2, y_intersect/2, 'VR');
    text((x_intersect+VS)/2, y_intersect/2, 'VZL');
    text(VS/2, -0.1, 'VS');
    text(x_intersect - y_intersect*current_length/2, y_intersect + x_intersect*current_length/2, 'I');

    % Title and labels
    title('Phasor Diagram with Circle Equations');
    xlabel('Real Component');
    ylabel('Imaginary Component');

    % Add origin point
    plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
    text(0.1, 0.1, 'O');

    % Add VS point
    plot(VS, 0, 'ko', 'MarkerFaceColor', 'k');
    text(VS+0.1, 0.1, 'VS');

    % Add legend
    legend('VR Circle', 'VZL Circle', 'VR Fasor', 'VZL Fasor', 'VS Fasor', 'Current I', 'Intersection');
end

% Function to calculate and display impedance
function calculate_and_display_impedance(VS, VR, VZL, R, angle_VZL)
    % Calculate current I = VR/R
    I = VR/R;

    % Impedance ZL = VZL/I with the calculated angle
    ZL_magnitude = VZL / I;

    % For an RC series circuit in the positive quadrant
    ZL_real = ZL_magnitude * cos(angle_VZL);
    ZL_imag = -ZL_magnitude * sin(angle_VZL); % Negative for capacitive impedance

    % Complete impedance
    Z = ZL_real + ZL_imag * 1i;

    % RC components calculation
    R_component = ZL_real;

    % For a capacitive reactance Xc = 1/(ωC)
    % Assuming ω = 1 rad/s for simplicity
    if ZL_imag ~= 0
        C_component = -1/ZL_imag;
    else
        C_component = Inf;
    end

    % Display results
    fprintf('\nResults:\n');
    fprintf('-------\n');
    fprintf('Load Impedance (ZL) = %.4f + j(%.4f) Ohms\n', real(Z), imag(Z));
    fprintf('Magnitude of ZL = %.4f Ohms\n', abs(Z));
    fprintf('Phase angle = %.4f degrees\n', angle(Z) * 180/pi);
    fprintf('Equivalent RC components:\n');
    fprintf('  R = %.4f Ohms\n', R_component);
    fprintf('  C = %.4e Farads (assuming ω = 1 rad/s)\n', C_component);

    fprintf('\nNote: Use the zoom tools in the figure window to examine the phasor diagram\n');
end
